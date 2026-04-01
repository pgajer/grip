#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <thread>
#include <unordered_map>
#include <limits>
#include <queue>
#include <string>
#include <vector>

#include "Graph.h"
#include "Point.h"

namespace {

struct GeodesicPathEdge {
    size_tt u;
    size_tt v;
    double coeff;
};

struct GeodesicPairCache {
    coord_t graphDistance;
    std::vector<GeodesicPathEdge> pathEdges;
};

struct GeodesicMdsState {
    double energy;
    double gmdsEnergy;
    double anchorEnergy;
    double edgeSpringEnergy;
    double repulsionEnergy;
    int repulsionPairCount;
    int repulsionActivePairCount;
    double smoothEnergy;
    double gradNorm2;
    std::vector<Point<>> gradient;
};

struct FlatGeodesicCacheView {
    std::vector<int> pairOffsets;
    std::vector<int> edgeU;
    std::vector<int> edgeV;
    std::vector<double> edgeCoeff;
    std::vector<double> pairGraphDistance;

    size_t pairCount() const
    {
        return pairGraphDistance.size();
    }
};

struct FlatSmoothnessView {
    std::vector<int> offsets;
    std::vector<int> neighbors;

    bool enabled() const
    {
        return !offsets.empty();
    }

    size_t vertexCount() const
    {
        return offsets.empty() ? 0 : offsets.size() - 1;
    }
};

struct FlatEdgePenaltyView {
    std::vector<int> u;
    std::vector<int> v;
    std::vector<double> target;

    bool enabled() const
    {
        return !u.empty();
    }

    size_t edgeCount() const
    {
        return u.size();
    }
};

struct FlatRepulsionView {
    std::vector<int> u;
    std::vector<int> v;
    std::vector<double> target;

    bool enabled() const
    {
        return !u.empty();
    }

    size_t pairCount() const
    {
        return u.size();
    }
};

void validate_geodesic_mds_args(int max_iter,
                                double edge_length_epsilon,
                                double initial_step,
                                double step_shrink,
                                double armijo_factor,
                                double grad_tol,
                                double min_step)
{
    if(max_iter < 0)
        Rcpp::stop("max_iter must be a non-negative integer");
    if(!std::isfinite(edge_length_epsilon) || edge_length_epsilon < 0.0)
        Rcpp::stop("edge_length_epsilon must be finite and >= 0");
    if(!std::isfinite(initial_step) || initial_step <= 0.0)
        Rcpp::stop("initial_step must be finite and > 0");
    if(!std::isfinite(step_shrink) || step_shrink <= 0.0 || step_shrink >= 1.0)
        Rcpp::stop("step_shrink must be finite and in (0, 1)");
    if(!std::isfinite(armijo_factor) || armijo_factor < 0.0)
        Rcpp::stop("armijo_factor must be finite and >= 0");
    if(!std::isfinite(grad_tol) || grad_tol < 0.0)
        Rcpp::stop("grad_tol must be finite and >= 0");
    if(!std::isfinite(min_step) || min_step <= 0.0)
        Rcpp::stop("min_step must be finite and > 0");
}

int resolve_gmds_thread_count(int requested,
                              size_t pairCount)
{
    if(pairCount == 0)
        return 1;

    if(requested > 0)
        return std::max(1, requested);

    const char *envValue = std::getenv("GRIP_GMDS_THREADS");
    if(envValue != nullptr){
        int parsed = std::atoi(envValue);
        if(parsed > 0)
            return parsed;
    }

    unsigned int hw = std::thread::hardware_concurrency();
    if(hw == 0)
        return 1;
    return static_cast<int>(std::max<unsigned int>(1U, hw));
}

double logsumexp_values(const std::vector<double> &values)
{
    double maxValue = -std::numeric_limits<double>::infinity();
    for(size_t i = 0; i < values.size(); i++){
        if(values[i] > maxValue)
            maxValue = values[i];
    }
    if(!std::isfinite(maxValue))
        return maxValue;

    double total = 0.0;
    for(size_t i = 0; i < values.size(); i++)
        total += std::exp(values[i] - maxValue);
    return maxValue + std::log(total);
}

void populate_graph_from_adj_list(Graph &graph,
                                  Rcpp::List adj_list,
                                  Rcpp::Nullable<Rcpp::List> weight_list,
                                  int n)
{
    if(n <= 0)
        Rcpp::stop("n must be positive");
    if(adj_list.size() != n)
        Rcpp::stop("adj_list length must match n");

    bool useWeights = weight_list.isNotNull();
    Rcpp::List weight_list_val;
    if(useWeights){
        weight_list_val = weight_list.get();
        if(weight_list_val.size() != n)
            Rcpp::stop("weight_list length must match n");
    }

    std::vector<std::vector<size_tt>> adj(static_cast<size_t>(n));
    std::vector<std::vector<coord_t>> weights;
    if(useWeights)
        weights.resize(static_cast<size_t>(n));

    for(int i = 0; i < n; i++){
        Rcpp::IntegerVector neigh = adj_list[i];
        adj[static_cast<size_t>(i)].reserve(static_cast<size_t>(neigh.size()));
        Rcpp::NumericVector w;
        if(useWeights){
            w = weight_list_val[i];
            if(neigh.size() != w.size())
                Rcpp::stop("weight_list must be parallel to adj_list");
            weights[static_cast<size_t>(i)].reserve(static_cast<size_t>(w.size()));
        }
        for(int j = 0; j < neigh.size(); j++){
            int v = neigh[j];
            if(v <= 0 || v > n)
                Rcpp::stop("adj_list must be 1-based and within [1, n]");
            adj[static_cast<size_t>(i)].push_back(static_cast<size_tt>(v - 1));
            if(useWeights){
                double wj = w[j];
                if(!std::isfinite(wj) || wj <= 0.0)
                    Rcpp::stop("weight_list must contain finite values > 0; invalid value %.17g at weight_list[[%d]][%d]",
                               wj,
                               i + 1,
                               j + 1);
                weights[static_cast<size_t>(i)].push_back(static_cast<coord_t>(wj));
            }
        }
    }

    graph.from_adj_list(static_cast<size_tt>(n), adj, useWeights ? &weights : nullptr);
}

std::vector<Point<>> matrix_to_points(const Rcpp::NumericMatrix &coords)
{
    int n = coords.nrow();
    int dim = coords.ncol();
    std::vector<Point<>> out(static_cast<size_t>(n));
    for(int i = 0; i < n; i++){
        out[static_cast<size_t>(i)].set_to_zero();
        out[static_cast<size_t>(i)].setX(coords(i, 0));
        if(dim > 1)
            out[static_cast<size_t>(i)].setY(coords(i, 1));
        if(dim > 2)
            out[static_cast<size_t>(i)].setZ(coords(i, 2));
    }
    return out;
}

Rcpp::NumericMatrix points_to_matrix(const std::vector<Point<>> &coords,
                                     int dim)
{
    Rcpp::NumericMatrix out(static_cast<int>(coords.size()), dim);
    for(int i = 0; i < out.nrow(); i++){
        out(i, 0) = coords[static_cast<size_t>(i)].getX();
        if(dim > 1)
            out(i, 1) = coords[static_cast<size_t>(i)].getY();
        if(dim > 2)
            out(i, 2) = coords[static_cast<size_t>(i)].getZ();
    }
    return out;
}

void recenter_points(std::vector<Point<>> &coords)
{
    if(coords.empty())
        return;
    Point<> mean;
    mean.set_to_zero();
    for(size_t i = 0; i < coords.size(); i++)
        mean += coords[i];
    mean /= static_cast<double>(coords.size());
    for(size_t i = 0; i < coords.size(); i++)
        coords[i] -= mean;
}

void compute_shortest_path_tree(const Graph &graph,
                                size_tt source,
                                std::vector<double> &dist,
                                std::vector<int> *parent)
{
    size_tt n = graph.get_numOfVert();
    const double inf = std::numeric_limits<double>::infinity();
    const double tol = 1e-10;
    dist.assign(static_cast<size_t>(n), inf);
    if(parent)
        parent->assign(static_cast<size_t>(n), -1);
    if(source >= n)
        return;

    dist[static_cast<size_t>(source)] = 0.0;

    if(!graph.has_weights()){
        std::queue<size_tt> q;
        q.push(source);

        while(!q.empty()){
            size_tt current = q.front();
            q.pop();
            double currentDist = dist[static_cast<size_t>(current)];
            size_tt deg = graph.get_Deg(current);
            for(size_tt adjIndex = 0; adjIndex < deg; adjIndex++){
                size_tt overt = graph.get_adjVert(current, adjIndex);
                double alt = currentDist + 1.0;
                double best = dist[static_cast<size_t>(overt)];
                bool improve = alt + tol < best;
                bool equal = std::isfinite(best) && std::fabs(alt - best) <= tol;
                if(improve){
                    dist[static_cast<size_t>(overt)] = alt;
                    if(parent)
                        (*parent)[static_cast<size_t>(overt)] = static_cast<int>(current);
                    q.push(overt);
                } else if(equal && parent){
                    int currentParent = (*parent)[static_cast<size_t>(overt)];
                    if(currentParent < 0 || current < static_cast<size_tt>(currentParent))
                        (*parent)[static_cast<size_t>(overt)] = static_cast<int>(current);
                }
            }
        }
        return;
    }

    struct QueueNode {
        double dist;
        size_tt vert;
    };
    struct QueueNodeGreater {
        bool operator()(const QueueNode &lhs, const QueueNode &rhs) const
        {
            if(lhs.dist != rhs.dist)
                return lhs.dist > rhs.dist;
            return lhs.vert > rhs.vert;
        }
    };
    std::priority_queue<QueueNode,
                        std::vector<QueueNode>,
                        QueueNodeGreater> pq;
    pq.push(QueueNode{0.0, source});

    while(!pq.empty()){
        double currentDist = pq.top().dist;
        size_tt current = pq.top().vert;
        pq.pop();
        if(currentDist > dist[static_cast<size_t>(current)] + tol)
            continue;
        size_tt deg = graph.get_Deg(current);
        for(size_tt adjIndex = 0; adjIndex < deg; adjIndex++){
            size_tt overt = graph.get_adjVert(current, adjIndex);
            double alt = currentDist + graph.get_edge_weight(current, adjIndex);
            double best = dist[static_cast<size_t>(overt)];
            double scale = std::max(1.0,
                                    std::max(std::fabs(alt),
                                             std::isfinite(best) ? std::fabs(best) : 0.0));
            bool improve = !std::isfinite(best) || alt + tol * scale < best;
            bool equal = std::isfinite(best) &&
                std::fabs(alt - best) <= tol * scale;
            if(improve){
                dist[static_cast<size_t>(overt)] = alt;
                if(parent)
                    (*parent)[static_cast<size_t>(overt)] = static_cast<int>(current);
                pq.push(QueueNode{alt, overt});
            } else if(equal && parent){
                int currentParent = (*parent)[static_cast<size_t>(overt)];
                if(currentParent < 0 || current < static_cast<size_tt>(currentParent))
                    (*parent)[static_cast<size_t>(overt)] = static_cast<int>(current);
            }
        }
    }
}

std::vector<size_tt> ordered_vertices_by_distance(const std::vector<double> &distRow)
{
    std::vector<size_tt> order(distRow.size());
    std::iota(order.begin(), order.end(), static_cast<size_tt>(0));
    std::sort(order.begin(), order.end(), [&distRow](size_tt lhs, size_tt rhs) {
        double dl = distRow[static_cast<size_t>(lhs)];
        double dr = distRow[static_cast<size_t>(rhs)];
        bool lf = std::isfinite(dl);
        bool rf = std::isfinite(dr);
        if(lf != rf)
            return lf > rf;
        if(dl != dr)
            return dl < dr;
        return lhs < rhs;
    });
    return order;
}

void compute_shortest_path_predecessors(const Graph &graph,
                                        size_tt source,
                                        const std::vector<double> &distRow,
                                        std::vector<std::vector<size_tt>> &predecessors)
{
    size_tt n = graph.get_numOfVert();
    const double tol = std::sqrt(std::numeric_limits<double>::epsilon());
    predecessors.assign(static_cast<size_t>(n), std::vector<size_tt>());
    if(source >= n)
        return;

    for(size_tt u = 0; u < n; u++){
        double dU = distRow[static_cast<size_t>(u)];
        if(!std::isfinite(dU))
            continue;
        size_tt deg = graph.get_Deg(u);
        for(size_tt adjIndex = 0; adjIndex < deg; adjIndex++){
            size_tt v = graph.get_adjVert(u, adjIndex);
            double dV = distRow[static_cast<size_t>(v)];
            if(!std::isfinite(dV))
                continue;
            double step = graph.has_weights() ? graph.get_edge_weight(u, adjIndex) : 1.0;
            double scale = std::max(1.0,
                                    std::max(std::fabs(dU),
                                             std::max(std::fabs(dV), std::fabs(step))));
            if(dU + tol * scale < dV &&
               std::fabs((dU + step) - dV) <= tol * scale)
                predecessors[static_cast<size_t>(v)].push_back(u);
        }
    }
}

std::vector<std::vector<size_tt>> compute_shortest_path_successors(
    const std::vector<std::vector<size_tt>> &predecessors)
{
    std::vector<std::vector<size_tt>> successors(predecessors.size());
    for(size_t v = 0; v < predecessors.size(); v++){
        for(size_t idx = 0; idx < predecessors[v].size(); idx++){
            size_tt u = predecessors[v][idx];
            successors[static_cast<size_t>(u)].push_back(static_cast<size_tt>(v));
        }
    }
    return successors;
}

std::vector<double> shortest_path_log_counts_forward(
    const std::vector<std::vector<size_tt>> &predecessors,
    size_tt source,
    const std::vector<size_tt> &orderVertices)
{
    std::vector<double> out(predecessors.size(),
                            -std::numeric_limits<double>::infinity());
    if(source >= predecessors.size())
        return out;
    out[static_cast<size_t>(source)] = 0.0;

    for(size_t idx = 0; idx < orderVertices.size(); idx++){
        size_tt v = orderVertices[idx];
        if(v == source)
            continue;
        const std::vector<size_tt> &pv = predecessors[static_cast<size_t>(v)];
        if(pv.empty())
            continue;
        std::vector<double> vals;
        vals.reserve(pv.size());
        for(size_t j = 0; j < pv.size(); j++)
            vals.push_back(out[static_cast<size_t>(pv[j])]);
        out[static_cast<size_t>(v)] = logsumexp_values(vals);
    }
    return out;
}

std::vector<int> shortest_path_ancestor_mask(
    const std::vector<std::vector<size_tt>> &predecessors,
    size_tt target)
{
    std::vector<int> keep(predecessors.size(), 0);
    if(target >= predecessors.size())
        return keep;

    std::queue<size_tt> q;
    q.push(target);
    keep[static_cast<size_t>(target)] = 1;
    while(!q.empty()){
        size_tt v = q.front();
        q.pop();
        const std::vector<size_tt> &pv = predecessors[static_cast<size_t>(v)];
        for(size_t idx = 0; idx < pv.size(); idx++){
            size_tt u = pv[idx];
            if(!keep[static_cast<size_t>(u)]){
                keep[static_cast<size_t>(u)] = 1;
                q.push(u);
            }
        }
    }
    return keep;
}

std::vector<double> shortest_path_log_counts_backward(
    const std::vector<std::vector<size_tt>> &successors,
    size_tt target,
    const std::vector<int> &ancestorMask,
    const std::vector<size_tt> &orderVertices)
{
    std::vector<double> out(successors.size(),
                            -std::numeric_limits<double>::infinity());
    if(target >= successors.size())
        return out;
    out[static_cast<size_t>(target)] = 0.0;

    for(size_t revIdx = orderVertices.size(); revIdx > 0; revIdx--){
        size_tt u = orderVertices[revIdx - 1];
        if(!ancestorMask[static_cast<size_t>(u)] || u == target)
            continue;
        const std::vector<size_tt> &sv = successors[static_cast<size_t>(u)];
        std::vector<double> vals;
        vals.reserve(sv.size());
        for(size_t idx = 0; idx < sv.size(); idx++){
            size_tt v = sv[idx];
            if(ancestorMask[static_cast<size_t>(v)] &&
               std::isfinite(out[static_cast<size_t>(v)]))
                vals.push_back(out[static_cast<size_t>(v)]);
        }
        if(!vals.empty())
            out[static_cast<size_t>(u)] = logsumexp_values(vals);
    }
    return out;
}

std::vector<GeodesicPairCache> build_all_pairs_cache(const Graph &graph)
{
    size_tt n = graph.get_numOfVert();
    std::vector<GeodesicPairCache> pairs;
    if(n < 2)
        return pairs;

    std::vector<double> dist;
    std::vector<int> parent;
    for(size_tt source = 0; source < n; source++){
        compute_shortest_path_tree(graph, source, dist, &parent);
        for(size_tt target = source + 1; target < n; target++){
            if(!std::isfinite(dist[static_cast<size_t>(target)]))
                Rcpp::stop("geodesic-MDS currently requires a connected graph");

            std::vector<size_tt> pathVertices;
            int current = static_cast<int>(target);
            while(current >= 0 && static_cast<size_tt>(current) != source){
                pathVertices.push_back(static_cast<size_tt>(current));
                current = parent[static_cast<size_t>(current)];
            }
            if(current < 0)
                Rcpp::stop("could not reconstruct a deterministic shortest path");
            pathVertices.push_back(source);
            std::reverse(pathVertices.begin(), pathVertices.end());

            GeodesicPairCache pair;
            pair.graphDistance = static_cast<coord_t>(dist[static_cast<size_t>(target)]);
            pair.pathEdges.reserve(pathVertices.size() > 0 ? pathVertices.size() - 1 : 0);
            for(size_t edgeIndex = 1; edgeIndex < pathVertices.size(); edgeIndex++){
                pair.pathEdges.push_back(
                    GeodesicPathEdge{pathVertices[edgeIndex - 1], pathVertices[edgeIndex], 1.0}
                );
            }
            pairs.push_back(std::move(pair));
        }
    }

    return pairs;
}

std::vector<GeodesicPairCache> build_cache_from_lists(Rcpp::List path_edges,
                                                      Rcpp::Nullable<Rcpp::List> path_edge_weights,
                                                      Rcpp::NumericVector pair_graph_distance,
                                                      int n)
{
    if(path_edges.size() != pair_graph_distance.size())
        Rcpp::stop("path_edges and pair_graph_distance must have the same length");

    bool useWeights = path_edge_weights.isNotNull();
    Rcpp::List weightList;
    if(useWeights){
        weightList = path_edge_weights.get();
        if(weightList.size() != path_edges.size())
            Rcpp::stop("path_edge_weights must be parallel to path_edges");
    }

    std::vector<GeodesicPairCache> pairs(static_cast<size_t>(path_edges.size()));
    for(R_xlen_t i = 0; i < path_edges.size(); i++){
        double graphDistance = pair_graph_distance[i];
        if(!std::isfinite(graphDistance) || graphDistance <= 0.0)
            Rcpp::stop("pair_graph_distance must contain finite values > 0");

        Rcpp::IntegerMatrix edges = path_edges[i];
        if(edges.ncol() != 2)
            Rcpp::stop("each element of path_edges must be a two-column integer matrix");

        Rcpp::NumericVector coeffs;
        if(useWeights){
            coeffs = weightList[i];
            if(coeffs.size() != 0 && coeffs.size() != edges.nrow())
                Rcpp::stop("path_edge_weights must be parallel to each path_edges matrix");
        }

        GeodesicPairCache pair;
        pair.graphDistance = graphDistance;
        pair.pathEdges.reserve(static_cast<size_t>(edges.nrow()));
        for(int row = 0; row < edges.nrow(); row++){
            int u = edges(row, 0);
            int v = edges(row, 1);
            if(u <= 0 || u > n || v <= 0 || v > n)
                Rcpp::stop("path_edges must use 1-based vertex ids within [1, nrow(coords)]");
            double coeff = coeffs.size() == 0 ? 1.0 : coeffs[row];
            if(!std::isfinite(coeff) || coeff < 0.0)
                Rcpp::stop("path_edge_weights must contain finite values >= 0");
            pair.pathEdges.push_back(
                GeodesicPathEdge{
                    static_cast<size_tt>(u - 1),
                    static_cast<size_tt>(v - 1),
                    coeff
                }
            );
        }
        pairs[static_cast<size_t>(i)] = std::move(pair);
    }

    return pairs;
}

Rcpp::List build_tie_average_shortest_path_cache_cpp_impl(
    Rcpp::List adj_list,
    Rcpp::Nullable<Rcpp::List> weight_list,
    Rcpp::IntegerMatrix pair_matrix,
    Rcpp::NumericMatrix dist_matrix)
{
    if(pair_matrix.ncol() != 2)
        Rcpp::stop("pair_matrix must have two columns");
    if(dist_matrix.nrow() != dist_matrix.ncol())
        Rcpp::stop("dist_matrix must be square");

    int n = dist_matrix.nrow();
    Graph graph;
    populate_graph_from_adj_list(graph, adj_list, weight_list, n);
    if(static_cast<int>(graph.get_numOfVert()) != n)
        Rcpp::stop("adj_list size must match nrow(dist_matrix)");

    int nPairs = pair_matrix.nrow();
    std::vector<std::vector<int>> pairIndicesBySource(static_cast<size_t>(n));
    for(int i = 0; i < nPairs; i++){
        int source = pair_matrix(i, 0);
        int target = pair_matrix(i, 1);
        if(source <= 0 || source > n || target <= 0 || target > n || source >= target)
            Rcpp::stop("pair_matrix must use 1-based ids with source < target");
        pairIndicesBySource[static_cast<size_t>(source - 1)].push_back(i);
    }

    std::vector<std::vector<int>> pairEdgeU(static_cast<size_t>(nPairs));
    std::vector<std::vector<int>> pairEdgeV(static_cast<size_t>(nPairs));
    std::vector<std::vector<double>> pairEdgeCoeff(static_cast<size_t>(nPairs));
    std::vector<double> pairGraphDistance(static_cast<size_t>(nPairs), NA_REAL);
    std::vector<double> pairPathCountLog(static_cast<size_t>(nPairs), NA_REAL);

    std::vector<double> distRow(static_cast<size_t>(n));
    std::vector<std::vector<size_tt>> predecessors;
    for(int source = 0; source < n; source++){
        if(pairIndicesBySource[static_cast<size_t>(source)].empty())
            continue;

        for(int col = 0; col < n; col++)
            distRow[static_cast<size_t>(col)] = dist_matrix(source, col);
        std::vector<size_tt> orderVertices = ordered_vertices_by_distance(distRow);
        compute_shortest_path_predecessors(graph,
                                           static_cast<size_tt>(source),
                                           distRow,
                                           predecessors);
        std::vector<std::vector<size_tt>> successors =
            compute_shortest_path_successors(predecessors);
        std::vector<double> logCountFrom =
            shortest_path_log_counts_forward(predecessors,
                                             static_cast<size_tt>(source),
                                             orderVertices);

        const std::vector<int> &pairIndices = pairIndicesBySource[static_cast<size_t>(source)];
        for(size_t idxPos = 0; idxPos < pairIndices.size(); idxPos++){
            int pairIndex = pairIndices[idxPos];
            int target = pair_matrix(pairIndex, 1) - 1;
            if(!std::isfinite(distRow[static_cast<size_t>(target)]))
                Rcpp::stop("geodesic-MDS currently requires a connected graph");

            std::vector<int> ancestorMask =
                shortest_path_ancestor_mask(predecessors, static_cast<size_tt>(target));
            std::vector<double> logCountTo =
                shortest_path_log_counts_backward(successors,
                                                 static_cast<size_tt>(target),
                                                 ancestorMask,
                                                 orderVertices);
            double totalLogCount = logCountFrom[static_cast<size_t>(target)];
            if(!std::isfinite(totalLogCount))
                Rcpp::stop("failed to count tied shortest paths");

            std::vector<int> edgeU;
            std::vector<int> edgeV;
            std::vector<double> edgeCoeff;
            for(size_t ordIdx = 0; ordIdx < orderVertices.size(); ordIdx++){
                size_tt u = orderVertices[ordIdx];
                if(!ancestorMask[static_cast<size_t>(u)])
                    continue;
                const std::vector<size_tt> &sv = successors[static_cast<size_t>(u)];
                for(size_t succIdx = 0; succIdx < sv.size(); succIdx++){
                    size_tt v = sv[succIdx];
                    if(!ancestorMask[static_cast<size_t>(v)] ||
                       !std::isfinite(logCountTo[static_cast<size_t>(v)]))
                        continue;
                    double coeff = std::exp(logCountFrom[static_cast<size_t>(u)] +
                                            logCountTo[static_cast<size_t>(v)] -
                                            totalLogCount);
                    edgeU.push_back(static_cast<int>(u) + 1);
                    edgeV.push_back(static_cast<int>(v) + 1);
                    edgeCoeff.push_back(coeff);
                }
            }

            pairEdgeU[static_cast<size_t>(pairIndex)] = std::move(edgeU);
            pairEdgeV[static_cast<size_t>(pairIndex)] = std::move(edgeV);
            pairEdgeCoeff[static_cast<size_t>(pairIndex)] = std::move(edgeCoeff);
            pairGraphDistance[static_cast<size_t>(pairIndex)] =
                distRow[static_cast<size_t>(target)];
            pairPathCountLog[static_cast<size_t>(pairIndex)] = totalLogCount;
        }
    }

    Rcpp::List pathEdges(nPairs);
    Rcpp::List pathEdgeWeights(nPairs);
    std::vector<int> flatPairOffsets;
    std::vector<int> flatEdgeU;
    std::vector<int> flatEdgeV;
    std::vector<double> flatEdgeCoeff;
    flatPairOffsets.reserve(static_cast<size_t>(nPairs + 1));
    flatPairOffsets.push_back(0);
    for(int pairIndex = 0; pairIndex < nPairs; pairIndex++){
        const std::vector<int> &edgeU = pairEdgeU[static_cast<size_t>(pairIndex)];
        const std::vector<int> &edgeV = pairEdgeV[static_cast<size_t>(pairIndex)];
        const std::vector<double> &edgeCoeff = pairEdgeCoeff[static_cast<size_t>(pairIndex)];
        Rcpp::IntegerMatrix edgeMat(static_cast<int>(edgeU.size()), 2);
        Rcpp::NumericVector coeffVec(static_cast<int>(edgeCoeff.size()));
        for(size_t idx = 0; idx < edgeU.size(); idx++){
            edgeMat(static_cast<int>(idx), 0) = edgeU[idx];
            edgeMat(static_cast<int>(idx), 1) = edgeV[idx];
            coeffVec[static_cast<int>(idx)] = edgeCoeff[idx];
            flatEdgeU.push_back(edgeU[idx] - 1);
            flatEdgeV.push_back(edgeV[idx] - 1);
            flatEdgeCoeff.push_back(edgeCoeff[idx]);
        }
        flatPairOffsets.push_back(static_cast<int>(flatEdgeU.size()));
        pathEdges[pairIndex] = edgeMat;
        pathEdgeWeights[pairIndex] = coeffVec;
    }

    return Rcpp::List::create(
        Rcpp::_["path_edges"] = pathEdges,
        Rcpp::_["path_edge_weights"] = pathEdgeWeights,
        Rcpp::_["pair_graph_distance"] = Rcpp::wrap(pairGraphDistance),
        Rcpp::_["pair_path_count_log"] = Rcpp::wrap(pairPathCountLog),
        Rcpp::_["flat_pair_edge_offsets"] = Rcpp::wrap(flatPairOffsets),
        Rcpp::_["flat_edge_u"] = Rcpp::wrap(flatEdgeU),
        Rcpp::_["flat_edge_v"] = Rcpp::wrap(flatEdgeV),
        Rcpp::_["flat_edge_coeff"] = Rcpp::wrap(flatEdgeCoeff)
    );
}

std::vector<double> resolve_anchor_schedule(Rcpp::Nullable<Rcpp::NumericVector> anchor_weights,
                                            int max_iter)
{
    std::vector<double> out(static_cast<size_t>(max_iter + 1), 0.0);
    if(anchor_weights.isNull())
        return out;

    Rcpp::NumericVector raw = anchor_weights.get();
    if(raw.size() == 0)
        return out;
    if(raw.size() != 1 && raw.size() != max_iter + 1)
        Rcpp::stop("anchor_weights must have length 1 or max_iter + 1");

    for(int i = 0; i <= max_iter; i++){
        double value = raw.size() == 1 ? raw[0] : raw[i];
        if(!std::isfinite(value) || value < 0.0)
            Rcpp::stop("anchor_weights must contain finite values >= 0");
        out[static_cast<size_t>(i)] = value;
    }
    return out;
}

FlatSmoothnessView build_flat_smoothness_view(Rcpp::IntegerVector smooth_adj_offsets,
                                              Rcpp::IntegerVector smooth_adj_vertices,
                                              int n)
{
    FlatSmoothnessView view{
        Rcpp::as<std::vector<int>>(smooth_adj_offsets),
        Rcpp::as<std::vector<int>>(smooth_adj_vertices)
    };
    if(view.offsets.empty() && view.neighbors.empty())
        return view;

    if(static_cast<int>(view.offsets.size()) != n + 1)
        Rcpp::stop("smooth_adj_offsets must have length nrow(coords) + 1");
    if(view.offsets.front() != 0)
        Rcpp::stop("smooth_adj_offsets must start at 0");
    if(view.offsets.back() != static_cast<int>(view.neighbors.size()))
        Rcpp::stop("smooth_adj_offsets must end at length(smooth_adj_vertices)");
    for(int i = 0; i < n; i++){
        if(view.offsets[static_cast<size_t>(i)] > view.offsets[static_cast<size_t>(i + 1)])
            Rcpp::stop("smooth_adj_offsets must be nondecreasing");
    }
    for(size_t i = 0; i < view.neighbors.size(); i++){
        int v = view.neighbors[i];
        if(v < 0 || v >= n)
            Rcpp::stop("smooth_adj_vertices must contain 0-based vertex ids within [0, nrow(coords) - 1]");
    }
    return view;
}

GeodesicMdsState evaluate_state(const std::vector<Point<>> &coords,
                                const std::vector<GeodesicPairCache> &pairs,
                                double eps2,
                                const std::vector<Point<>> *anchor,
                                double anchorWeight)
{
    GeodesicMdsState state;
    state.energy = 0.0;
    state.gmdsEnergy = 0.0;
    state.anchorEnergy = 0.0;
    state.edgeSpringEnergy = 0.0;
    state.repulsionEnergy = 0.0;
    state.smoothEnergy = 0.0;
    state.gradNorm2 = 0.0;
    state.gradient.assign(coords.size(), Point<>());
    for(size_t i = 0; i < state.gradient.size(); i++)
        state.gradient[i].set_to_zero();

    for(size_t pairIndex = 0; pairIndex < pairs.size(); pairIndex++){
        const GeodesicPairCache &pair = pairs[pairIndex];
        std::vector<Point<>> edgeDiffs;
        std::vector<double> edgeLens;
        edgeDiffs.reserve(pair.pathEdges.size());
        edgeLens.reserve(pair.pathEdges.size());

        double h = 0.0;
        for(size_t edgeIndex = 0; edgeIndex < pair.pathEdges.size(); edgeIndex++){
            const GeodesicPathEdge &edgeRef = pair.pathEdges[edgeIndex];
            Point<> diff = coords[static_cast<size_t>(edgeRef.u)] -
                           coords[static_cast<size_t>(edgeRef.v)];
            double len = std::sqrt(diff.fnorm2() + eps2);
            edgeDiffs.push_back(diff);
            edgeLens.push_back(len);
            h += edgeRef.coeff * len;
        }

        double resid = h - pair.graphDistance;
        state.gmdsEnergy += 0.5 * resid * resid;
        for(size_t edgeIndex = 0; edgeIndex < pair.pathEdges.size(); edgeIndex++){
            if(edgeLens[edgeIndex] <= 0.0)
                continue;
            const GeodesicPathEdge &edgeRef = pair.pathEdges[edgeIndex];
            Point<> stepVec = edgeDiffs[edgeIndex] * (resid * edgeRef.coeff / edgeLens[edgeIndex]);
            state.gradient[static_cast<size_t>(edgeRef.u)] += stepVec;
            state.gradient[static_cast<size_t>(edgeRef.v)] -= stepVec;
        }
    }

    if(anchor && anchorWeight > 0.0){
        double rawPenalty = 0.0;
        for(size_t i = 0; i < coords.size(); i++){
            Point<> diff = coords[i] - (*anchor)[i];
            rawPenalty += diff.fnorm2();
            state.gradient[i] += diff * (2.0 * anchorWeight);
        }
        state.anchorEnergy = anchorWeight * rawPenalty;
    }

    state.energy = state.gmdsEnergy + state.anchorEnergy;

    for(size_t i = 0; i < state.gradient.size(); i++)
        state.gradNorm2 += state.gradient[i].fnorm2();

    return state;
}

void accumulate_flat_pair_range(const std::vector<Point<>> &coords,
                                const FlatGeodesicCacheView &cache,
                                double eps2,
                                size_t pairBegin,
                                size_t pairEnd,
                                std::vector<Point<>> &gradient,
                                double &gmdsEnergy)
{
    gmdsEnergy = 0.0;
    for(size_t i = 0; i < gradient.size(); i++)
        gradient[i].set_to_zero();

    for(size_t pairIndex = pairBegin; pairIndex < pairEnd; pairIndex++){
        int edgeBegin = cache.pairOffsets[pairIndex];
        int edgeEnd = cache.pairOffsets[pairIndex + 1];
        double h = 0.0;
        for(int edgeIndex = edgeBegin; edgeIndex < edgeEnd; edgeIndex++){
            int u = cache.edgeU[edgeIndex];
            int v = cache.edgeV[edgeIndex];
            Point<> diff = coords[static_cast<size_t>(u)] -
                           coords[static_cast<size_t>(v)];
            double len = std::sqrt(diff.fnorm2() + eps2);
            h += cache.edgeCoeff[edgeIndex] * len;
        }

        double resid = h - cache.pairGraphDistance[pairIndex];
        gmdsEnergy += 0.5 * resid * resid;
        for(int edgeIndex = edgeBegin; edgeIndex < edgeEnd; edgeIndex++){
            int u = cache.edgeU[edgeIndex];
            int v = cache.edgeV[edgeIndex];
            Point<> diff = coords[static_cast<size_t>(u)] -
                           coords[static_cast<size_t>(v)];
            double len = std::sqrt(diff.fnorm2() + eps2);
            if(len <= 0.0)
                continue;
            Point<> stepVec = diff * (resid * cache.edgeCoeff[edgeIndex] / len);
            gradient[static_cast<size_t>(u)] += stepVec;
            gradient[static_cast<size_t>(v)] -= stepVec;
        }
    }
}

void accumulate_flat_smoothness(const std::vector<Point<>> &coords,
                                const FlatSmoothnessView &smoothness,
                                double smoothWeight,
                                std::vector<Point<>> &gradient,
                                double &smoothEnergy)
{
    smoothEnergy = 0.0;
    if(!smoothness.enabled() || smoothWeight <= 0.0)
        return;

    std::vector<Point<>> residual(coords.size());
    for(size_t i = 0; i < residual.size(); i++)
        residual[i].set_to_zero();

    for(size_t i = 0; i < coords.size(); i++){
        int begin = smoothness.offsets[i];
        int end = smoothness.offsets[i + 1];
        int degree = end - begin;
        if(degree <= 0)
            continue;

        Point<> average;
        average.set_to_zero();
        for(int edgeIndex = begin; edgeIndex < end; edgeIndex++)
            average += coords[static_cast<size_t>(smoothness.neighbors[edgeIndex])];
        average /= static_cast<double>(degree);

        residual[i] = coords[i] - average;
        smoothEnergy += residual[i].fnorm2();
        gradient[i] += residual[i] * (2.0 * smoothWeight);
    }

    for(size_t i = 0; i < coords.size(); i++){
        int begin = smoothness.offsets[i];
        int end = smoothness.offsets[i + 1];
        int degree = end - begin;
        if(degree <= 0)
            continue;

        Point<> shared = residual[i] * (2.0 * smoothWeight / static_cast<double>(degree));
        for(int edgeIndex = begin; edgeIndex < end; edgeIndex++)
            gradient[static_cast<size_t>(smoothness.neighbors[edgeIndex])] -= shared;
    }

    smoothEnergy *= smoothWeight;
}

void accumulate_flat_edge_springs(const std::vector<Point<>> &coords,
                                  const FlatEdgePenaltyView &edgePenalty,
                                  double eps2,
                                  double edgeSpringWeight,
                                  std::vector<Point<>> &gradient,
                                  double &edgeSpringEnergy)
{
    edgeSpringEnergy = 0.0;
    if(!edgePenalty.enabled() || edgeSpringWeight <= 0.0)
        return;

    double rawPenalty = 0.0;
    for(size_t edgeIndex = 0; edgeIndex < edgePenalty.edgeCount(); edgeIndex++){
        int u = edgePenalty.u[edgeIndex];
        int v = edgePenalty.v[edgeIndex];
        Point<> diff = coords[static_cast<size_t>(u)] -
                       coords[static_cast<size_t>(v)];
        double len = std::sqrt(diff.fnorm2() + eps2);
        double resid = len - edgePenalty.target[edgeIndex];
        rawPenalty += resid * resid;
        if(len <= 0.0)
            continue;
        Point<> stepVec = diff * (edgeSpringWeight * resid / len);
        gradient[static_cast<size_t>(u)] += stepVec;
        gradient[static_cast<size_t>(v)] -= stepVec;
    }

    edgeSpringEnergy = 0.5 * edgeSpringWeight * rawPenalty;
}

void accumulate_flat_repulsion(const std::vector<Point<>> &coords,
                               const FlatRepulsionView &repulsion,
                               double eps2,
                               double repulsionWeight,
                               std::vector<Point<>> &gradient,
                               double &repulsionEnergy,
                               int &activePairCount)
{
    repulsionEnergy = 0.0;
    activePairCount = 0;
    if(!repulsion.enabled() || repulsionWeight <= 0.0)
        return;

    double rawPenalty = 0.0;
    for(size_t pairIndex = 0; pairIndex < repulsion.pairCount(); pairIndex++){
        int u = repulsion.u[pairIndex];
        int v = repulsion.v[pairIndex];
        Point<> diff = coords[static_cast<size_t>(u)] -
                       coords[static_cast<size_t>(v)];
        double len = std::sqrt(diff.fnorm2() + eps2);
        double resid = repulsion.target[pairIndex] - len;
        if(resid <= 0.0)
            continue;
        activePairCount++;
        rawPenalty += resid * resid;
        if(len <= 0.0)
            continue;
        Point<> stepVec = diff * (-repulsionWeight * resid / len);
        gradient[static_cast<size_t>(u)] += stepVec;
        gradient[static_cast<size_t>(v)] -= stepVec;
    }

    repulsionEnergy = 0.5 * repulsionWeight * rawPenalty;
}

GeodesicMdsState evaluate_flat_state(const std::vector<Point<>> &coords,
                                     const FlatGeodesicCacheView &cache,
                                     const FlatSmoothnessView &smoothness,
                                     const FlatEdgePenaltyView &edgePenalty,
                                     const FlatRepulsionView &repulsion,
                                     double eps2,
                                     const std::vector<Point<>> *anchor,
                                     double anchorWeight,
                                     double smoothWeight,
                                     double edgeSpringWeight,
                                     double repulsionWeight,
                                     int requestedThreads)
{
    GeodesicMdsState state;
    state.energy = 0.0;
    state.gmdsEnergy = 0.0;
    state.anchorEnergy = 0.0;
    state.edgeSpringEnergy = 0.0;
    state.repulsionEnergy = 0.0;
    state.repulsionPairCount = static_cast<int>(repulsion.pairCount());
    state.repulsionActivePairCount = 0;
    state.smoothEnergy = 0.0;
    state.gradNorm2 = 0.0;
    state.gradient.assign(coords.size(), Point<>());
    for(size_t i = 0; i < state.gradient.size(); i++)
        state.gradient[i].set_to_zero();

    int nThreads = resolve_gmds_thread_count(requestedThreads, cache.pairCount());
    if(nThreads > static_cast<int>(cache.pairCount()))
        nThreads = static_cast<int>(cache.pairCount());
    if(nThreads < 1)
        nThreads = 1;
    if(cache.pairCount() < 2048)
        nThreads = 1;

    if(nThreads == 1){
        accumulate_flat_pair_range(coords,
                                   cache,
                                   eps2,
                                   0,
                                   cache.pairCount(),
                                   state.gradient,
                                   state.gmdsEnergy);
    } else {
        std::vector<std::vector<Point<>>> threadGradients(
            static_cast<size_t>(nThreads),
            std::vector<Point<>>(coords.size())
        );
        std::vector<double> threadEnergy(static_cast<size_t>(nThreads), 0.0);
        std::vector<std::thread> workers;
        workers.reserve(static_cast<size_t>(nThreads));
        size_t pairCount = cache.pairCount();
        size_t chunk = (pairCount + static_cast<size_t>(nThreads) - 1) /
            static_cast<size_t>(nThreads);

        for(int threadIndex = 0; threadIndex < nThreads; threadIndex++){
            size_t begin = static_cast<size_t>(threadIndex) * chunk;
            size_t end = std::min(pairCount, begin + chunk);
            if(begin >= end)
                continue;
            workers.emplace_back([&, threadIndex, begin, end]() {
                accumulate_flat_pair_range(coords,
                                           cache,
                                           eps2,
                                           begin,
                                           end,
                                           threadGradients[static_cast<size_t>(threadIndex)],
                                           threadEnergy[static_cast<size_t>(threadIndex)]);
            });
        }
        for(size_t idx = 0; idx < workers.size(); idx++)
            workers[idx].join();

        for(int threadIndex = 0; threadIndex < nThreads; threadIndex++){
            state.gmdsEnergy += threadEnergy[static_cast<size_t>(threadIndex)];
            const std::vector<Point<>> &localGrad =
                threadGradients[static_cast<size_t>(threadIndex)];
            for(size_t vert = 0; vert < state.gradient.size(); vert++)
                state.gradient[vert] += localGrad[vert];
        }
    }

    accumulate_flat_smoothness(coords,
                               smoothness,
                               smoothWeight,
                               state.gradient,
                               state.smoothEnergy);
    accumulate_flat_edge_springs(coords,
                                 edgePenalty,
                                 eps2,
                                 edgeSpringWeight,
                                 state.gradient,
                                 state.edgeSpringEnergy);
    accumulate_flat_repulsion(coords,
                              repulsion,
                              eps2,
                              repulsionWeight,
                              state.gradient,
                              state.repulsionEnergy,
                              state.repulsionActivePairCount);

    if(anchor && anchorWeight > 0.0){
        double rawPenalty = 0.0;
        for(size_t i = 0; i < coords.size(); i++){
            Point<> diff = coords[i] - (*anchor)[i];
            rawPenalty += diff.fnorm2();
            state.gradient[i] += diff * (2.0 * anchorWeight);
        }
        state.anchorEnergy = anchorWeight * rawPenalty;
    }

    state.energy = state.gmdsEnergy + state.anchorEnergy + state.edgeSpringEnergy +
        state.repulsionEnergy + state.smoothEnergy;
    for(size_t i = 0; i < state.gradient.size(); i++)
        state.gradNorm2 += state.gradient[i].fnorm2();

    return state;
}

Rcpp::DataFrame build_trace_df(const std::vector<int> &iteration,
                               const std::vector<double> &energy,
                               const std::vector<double> &gmds_energy,
                               const std::vector<double> &anchor_energy,
                               const std::vector<double> &edge_spring_energy,
                               const std::vector<double> &repulsion_energy,
                               const std::vector<int> &repulsion_pair_count,
                               const std::vector<int> &repulsion_active_pair_count,
                               const std::vector<double> &smooth_energy,
                               const std::vector<double> &gradient_norm,
                               const std::vector<double> &step,
                               const std::vector<bool> &accepted,
                               const std::vector<double> &anchor_weight,
                               const std::vector<double> &edge_spring_weight,
                               const std::vector<double> &repulsion_weight,
                               const std::vector<double> &smooth_weight)
{
    return Rcpp::DataFrame::create(
        Rcpp::_["iteration"] = iteration,
        Rcpp::_["energy"] = energy,
        Rcpp::_["gmds_energy"] = gmds_energy,
        Rcpp::_["anchor_energy"] = anchor_energy,
        Rcpp::_["edge_spring_energy"] = edge_spring_energy,
        Rcpp::_["repulsion_energy"] = repulsion_energy,
        Rcpp::_["repulsion_pair_count"] = repulsion_pair_count,
        Rcpp::_["repulsion_active_pair_count"] = repulsion_active_pair_count,
        Rcpp::_["smooth_energy"] = smooth_energy,
        Rcpp::_["gradient_norm"] = gradient_norm,
        Rcpp::_["step"] = step,
        Rcpp::_["accepted"] = accepted,
        Rcpp::_["anchor_weight"] = anchor_weight,
        Rcpp::_["edge_spring_weight"] = edge_spring_weight,
        Rcpp::_["repulsion_weight"] = repulsion_weight,
        Rcpp::_["smooth_weight"] = smooth_weight,
        Rcpp::_["stringsAsFactors"] = false
    );
}

Rcpp::List build_frame_list(const std::vector<std::vector<Point<>>> &frames,
                            int dim)
{
    Rcpp::List out(frames.size());
    for(size_t i = 0; i < frames.size(); i++)
        out[static_cast<R_xlen_t>(i)] = points_to_matrix(frames[i], dim);
    return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List grip_build_tie_average_shortest_path_cache_cpp(
    Rcpp::List adj_list,
    Rcpp::Nullable<Rcpp::List> weight_list,
    Rcpp::IntegerMatrix pair_matrix,
    Rcpp::NumericMatrix dist_matrix)
{
    return build_tie_average_shortest_path_cache_cpp_impl(adj_list,
                                                          weight_list,
                                                          pair_matrix,
                                                          dist_matrix);
}

// [[Rcpp::export]]
Rcpp::List grip_optimize_geodesic_mds_adj_cpp(
    Rcpp::List adj_list,
    Rcpp::Nullable<Rcpp::List> weight_list,
    Rcpp::NumericMatrix coords,
    int max_iter,
    double edge_length_epsilon,
    double initial_step,
    double step_shrink,
    double armijo_factor,
    double grad_tol,
    double min_step,
    bool recenter,
    bool return_trace)
{
    if(coords.ncol() != 2 && coords.ncol() != 3)
        Rcpp::stop("coords must have 2 or 3 columns");
    if(coords.nrow() <= 1)
        Rcpp::stop("coords must have at least two rows");
    for(int i = 0; i < coords.size(); i++){
        if(!std::isfinite(coords[i]))
            Rcpp::stop("coords must contain only finite values");
    }

    validate_geodesic_mds_args(max_iter,
                               edge_length_epsilon,
                               initial_step,
                               step_shrink,
                               armijo_factor,
                               grad_tol,
                               min_step);

    Graph graph;
    populate_graph_from_adj_list(graph, adj_list, weight_list, coords.nrow());
    if(static_cast<int>(graph.get_numOfVert()) != coords.nrow())
        Rcpp::stop("adj_list size must match nrow(coords)");

    std::vector<GeodesicPairCache> pairs = build_all_pairs_cache(graph);
    std::vector<Point<>> current = matrix_to_points(coords);

    std::vector<int> trace_iteration;
    std::vector<double> trace_energy;
    std::vector<double> trace_gmds_energy;
    std::vector<double> trace_anchor_energy;
    std::vector<double> trace_edge_spring_energy;
    std::vector<double> trace_repulsion_energy;
    std::vector<int> trace_repulsion_pair_count;
    std::vector<int> trace_repulsion_active_pair_count;
    std::vector<double> trace_smooth_energy;
    std::vector<double> trace_gradient_norm;
    std::vector<double> trace_step;
    std::vector<bool> trace_accepted;
    std::vector<double> trace_anchor_weight;
    std::vector<double> trace_edge_spring_weight;
    std::vector<double> trace_repulsion_weight;
    std::vector<double> trace_smooth_weight;
    std::vector<std::vector<Point<>>> accepted_frames;
    accepted_frames.push_back(current);

    const double eps2 = edge_length_epsilon * edge_length_epsilon;
    const double gradTol2 = grad_tol * grad_tol;

    GeodesicMdsState state = evaluate_state(current, pairs, eps2, nullptr, 0.0);
    trace_iteration.push_back(0);
    trace_energy.push_back(state.energy);
    trace_gmds_energy.push_back(state.gmdsEnergy);
    trace_anchor_energy.push_back(state.anchorEnergy);
    trace_edge_spring_energy.push_back(0.0);
    trace_repulsion_energy.push_back(0.0);
    trace_repulsion_pair_count.push_back(0);
    trace_repulsion_active_pair_count.push_back(0);
    trace_smooth_energy.push_back(0.0);
    trace_gradient_norm.push_back(std::sqrt(state.gradNorm2));
    trace_step.push_back(NA_REAL);
    trace_accepted.push_back(true);
    trace_anchor_weight.push_back(0.0);
    trace_edge_spring_weight.push_back(0.0);
    trace_repulsion_weight.push_back(0.0);
    trace_smooth_weight.push_back(0.0);

    for(int iter = 1; iter <= max_iter; iter++){
        if(!std::isfinite(state.energy) || state.gradNorm2 <= gradTol2)
            break;

        double step = initial_step;
        bool accepted = false;
        std::vector<Point<>> proposal = current;
        GeodesicMdsState candidate = state;

        while(std::isfinite(step) && step >= min_step){
            proposal = current;
            for(size_t i = 0; i < proposal.size(); i++)
                proposal[i] -= state.gradient[i] * step;
            if(recenter)
                recenter_points(proposal);

            candidate = evaluate_state(proposal, pairs, eps2, nullptr, 0.0);
            double targetEnergy = state.energy - armijo_factor * step * state.gradNorm2;
            if(std::isfinite(candidate.energy) && candidate.energy <= targetEnergy){
                accepted = true;
                break;
            }
            step *= step_shrink;
        }

        trace_iteration.push_back(iter);
        trace_energy.push_back(accepted ? candidate.energy : state.energy);
        trace_gmds_energy.push_back(accepted ? candidate.gmdsEnergy : state.gmdsEnergy);
        trace_anchor_energy.push_back(accepted ? candidate.anchorEnergy : state.anchorEnergy);
        trace_edge_spring_energy.push_back(0.0);
        trace_repulsion_energy.push_back(0.0);
        trace_repulsion_pair_count.push_back(0);
        trace_repulsion_active_pair_count.push_back(0);
        trace_smooth_energy.push_back(0.0);
        trace_gradient_norm.push_back(std::sqrt(accepted ? candidate.gradNorm2 : state.gradNorm2));
        trace_step.push_back(accepted ? step : NA_REAL);
        trace_accepted.push_back(accepted);
        trace_anchor_weight.push_back(0.0);
        trace_edge_spring_weight.push_back(0.0);
        trace_repulsion_weight.push_back(0.0);
        trace_smooth_weight.push_back(0.0);

        if(!accepted)
            break;

        current = proposal;
        state = candidate;
        if(return_trace)
            accepted_frames.push_back(current);
    }

    if(!return_trace){
        accepted_frames.clear();
        accepted_frames.push_back(current);
    }

    return Rcpp::List::create(
        Rcpp::_["coords"] = points_to_matrix(current, coords.ncol()),
        Rcpp::_["trace"] = build_trace_df(trace_iteration,
                                          trace_energy,
                                          trace_gmds_energy,
                                          trace_anchor_energy,
                                          trace_edge_spring_energy,
                                          trace_repulsion_energy,
                                          trace_repulsion_pair_count,
                                          trace_repulsion_active_pair_count,
                                          trace_smooth_energy,
                                          trace_gradient_norm,
                                          trace_step,
                                          trace_accepted,
                                          trace_anchor_weight,
                                          trace_edge_spring_weight,
                                          trace_repulsion_weight,
                                          trace_smooth_weight),
        Rcpp::_["frames"] = build_frame_list(accepted_frames, coords.ncol())
    );
}

// [[Rcpp::export]]
Rcpp::List grip_optimize_geodesic_mds_cache_cpp(
    Rcpp::List path_edges,
    Rcpp::Nullable<Rcpp::List> path_edge_weights,
    Rcpp::NumericVector pair_graph_distance,
    Rcpp::NumericMatrix coords,
    int max_iter,
    double edge_length_epsilon,
    double initial_step,
    double step_shrink,
    double armijo_factor,
    double grad_tol,
    double min_step,
    bool recenter,
    bool return_trace,
    Rcpp::Nullable<Rcpp::NumericMatrix> anchor_coords = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> anchor_weights = R_NilValue)
{
    if(coords.ncol() != 2 && coords.ncol() != 3)
        Rcpp::stop("coords must have 2 or 3 columns");
    if(coords.nrow() <= 1)
        Rcpp::stop("coords must have at least two rows");
    for(int i = 0; i < coords.size(); i++){
        if(!std::isfinite(coords[i]))
            Rcpp::stop("coords must contain only finite values");
    }

    validate_geodesic_mds_args(max_iter,
                               edge_length_epsilon,
                               initial_step,
                               step_shrink,
                               armijo_factor,
                               grad_tol,
                               min_step);

    std::vector<GeodesicPairCache> pairs = build_cache_from_lists(
        path_edges,
        path_edge_weights,
        pair_graph_distance,
        coords.nrow()
    );
    std::vector<Point<>> current = matrix_to_points(coords);
    std::vector<Point<>> anchor;
    bool useAnchor = anchor_coords.isNotNull();
    if(useAnchor){
        Rcpp::NumericMatrix anchorMat = anchor_coords.get();
        if(anchorMat.nrow() != coords.nrow() || anchorMat.ncol() != coords.ncol())
            Rcpp::stop("anchor_coords must have the same dimensions as coords");
        for(int i = 0; i < anchorMat.size(); i++){
            if(!std::isfinite(anchorMat[i]))
                Rcpp::stop("anchor_coords must contain only finite values");
        }
        anchor = matrix_to_points(anchorMat);
    }
    std::vector<double> anchorSchedule = resolve_anchor_schedule(anchor_weights, max_iter);

    std::vector<int> trace_iteration;
    std::vector<double> trace_energy;
    std::vector<double> trace_gmds_energy;
    std::vector<double> trace_anchor_energy;
    std::vector<double> trace_edge_spring_energy;
    std::vector<double> trace_repulsion_energy;
    std::vector<int> trace_repulsion_pair_count;
    std::vector<int> trace_repulsion_active_pair_count;
    std::vector<double> trace_smooth_energy;
    std::vector<double> trace_gradient_norm;
    std::vector<double> trace_step;
    std::vector<bool> trace_accepted;
    std::vector<double> trace_anchor_weight;
    std::vector<double> trace_edge_spring_weight;
    std::vector<double> trace_repulsion_weight;
    std::vector<double> trace_smooth_weight;
    std::vector<std::vector<Point<>>> accepted_frames;
    accepted_frames.push_back(current);

    const double eps2 = edge_length_epsilon * edge_length_epsilon;
    const double gradTol2 = grad_tol * grad_tol;

    GeodesicMdsState state = evaluate_state(
        current,
        pairs,
        eps2,
        useAnchor ? &anchor : nullptr,
        anchorSchedule[0]
    );
    trace_iteration.push_back(0);
    trace_energy.push_back(state.energy);
    trace_gmds_energy.push_back(state.gmdsEnergy);
    trace_anchor_energy.push_back(state.anchorEnergy);
    trace_edge_spring_energy.push_back(0.0);
    trace_repulsion_energy.push_back(0.0);
    trace_repulsion_pair_count.push_back(0);
    trace_repulsion_active_pair_count.push_back(0);
    trace_smooth_energy.push_back(0.0);
    trace_gradient_norm.push_back(std::sqrt(state.gradNorm2));
    trace_step.push_back(NA_REAL);
    trace_accepted.push_back(true);
    trace_anchor_weight.push_back(anchorSchedule[0]);
    trace_edge_spring_weight.push_back(0.0);
    trace_repulsion_weight.push_back(0.0);
    trace_smooth_weight.push_back(0.0);

    for(int iter = 1; iter <= max_iter; iter++){
        double iterAnchorWeight = anchorSchedule[static_cast<size_t>(iter)];
        state = evaluate_state(
            current,
            pairs,
            eps2,
            useAnchor ? &anchor : nullptr,
            iterAnchorWeight
        );
        if(!std::isfinite(state.energy) || state.gradNorm2 <= gradTol2)
            break;

        double step = initial_step;
        bool accepted = false;
        std::vector<Point<>> proposal = current;
        GeodesicMdsState candidate = state;

        while(std::isfinite(step) && step >= min_step){
            proposal = current;
            for(size_t i = 0; i < proposal.size(); i++)
                proposal[i] -= state.gradient[i] * step;
            if(recenter)
                recenter_points(proposal);

            candidate = evaluate_state(
                proposal,
                pairs,
                eps2,
                useAnchor ? &anchor : nullptr,
                iterAnchorWeight
            );
            double targetEnergy = state.energy - armijo_factor * step * state.gradNorm2;
            if(std::isfinite(candidate.energy) && candidate.energy <= targetEnergy){
                accepted = true;
                break;
            }
            step *= step_shrink;
        }

        trace_iteration.push_back(iter);
        trace_energy.push_back(accepted ? candidate.energy : state.energy);
        trace_gmds_energy.push_back(accepted ? candidate.gmdsEnergy : state.gmdsEnergy);
        trace_anchor_energy.push_back(accepted ? candidate.anchorEnergy : state.anchorEnergy);
        trace_edge_spring_energy.push_back(0.0);
        trace_repulsion_energy.push_back(0.0);
        trace_repulsion_pair_count.push_back(0);
        trace_repulsion_active_pair_count.push_back(0);
        trace_smooth_energy.push_back(0.0);
        trace_gradient_norm.push_back(std::sqrt(accepted ? candidate.gradNorm2 : state.gradNorm2));
        trace_step.push_back(accepted ? step : NA_REAL);
        trace_accepted.push_back(accepted);
        trace_anchor_weight.push_back(iterAnchorWeight);
        trace_edge_spring_weight.push_back(0.0);
        trace_repulsion_weight.push_back(0.0);
        trace_smooth_weight.push_back(0.0);

        if(!accepted)
            break;

        current = proposal;
        state = candidate;
        if(return_trace)
            accepted_frames.push_back(current);
    }

    if(!return_trace){
        accepted_frames.clear();
        accepted_frames.push_back(current);
    }

    return Rcpp::List::create(
        Rcpp::_["coords"] = points_to_matrix(current, coords.ncol()),
        Rcpp::_["trace"] = build_trace_df(trace_iteration,
                                          trace_energy,
                                          trace_gmds_energy,
                                          trace_anchor_energy,
                                          trace_edge_spring_energy,
                                          trace_repulsion_energy,
                                          trace_repulsion_pair_count,
                                          trace_repulsion_active_pair_count,
                                          trace_smooth_energy,
                                          trace_gradient_norm,
                                          trace_step,
                                          trace_accepted,
                                          trace_anchor_weight,
                                          trace_edge_spring_weight,
                                          trace_repulsion_weight,
                                          trace_smooth_weight),
        Rcpp::_["frames"] = build_frame_list(accepted_frames, coords.ncol()),
        Rcpp::_["final_anchor_weight"] = trace_anchor_weight.back(),
        Rcpp::_["final_edge_spring_weight"] = trace_edge_spring_weight.back(),
        Rcpp::_["final_repulsion_weight"] = trace_repulsion_weight.back(),
        Rcpp::_["final_smoothness_weight"] = trace_smooth_weight.back()
    );
}

// [[Rcpp::export]]
Rcpp::List grip_optimize_geodesic_mds_flat_cpp(
    Rcpp::IntegerVector flat_pair_edge_offsets,
    Rcpp::IntegerVector flat_edge_u,
    Rcpp::IntegerVector flat_edge_v,
    Rcpp::NumericVector flat_edge_coeff,
    Rcpp::NumericVector pair_graph_distance,
    Rcpp::NumericMatrix coords,
    int max_iter,
    double edge_length_epsilon,
    double initial_step,
    double step_shrink,
    double armijo_factor,
    double grad_tol,
    double min_step,
    bool recenter,
    bool return_trace,
    Rcpp::Nullable<Rcpp::NumericMatrix> anchor_coords = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> anchor_weights = R_NilValue,
    Rcpp::IntegerVector smooth_adj_offsets = Rcpp::IntegerVector(),
    Rcpp::IntegerVector smooth_adj_vertices = Rcpp::IntegerVector(),
    Rcpp::Nullable<Rcpp::NumericVector> smooth_weights = R_NilValue,
    Rcpp::IntegerVector graph_edge_u = Rcpp::IntegerVector(),
    Rcpp::IntegerVector graph_edge_v = Rcpp::IntegerVector(),
    Rcpp::NumericVector graph_edge_target = Rcpp::NumericVector(),
    Rcpp::Nullable<Rcpp::NumericVector> edge_spring_weights = R_NilValue,
    Rcpp::IntegerVector repulsion_u = Rcpp::IntegerVector(),
    Rcpp::IntegerVector repulsion_v = Rcpp::IntegerVector(),
    Rcpp::NumericVector repulsion_target = Rcpp::NumericVector(),
    Rcpp::Nullable<Rcpp::NumericVector> repulsion_weights = R_NilValue,
    int n_threads = 0)
{
    if(coords.ncol() != 2 && coords.ncol() != 3)
        Rcpp::stop("coords must have 2 or 3 columns");
    if(coords.nrow() <= 1)
        Rcpp::stop("coords must have at least two rows");
    for(int i = 0; i < coords.size(); i++){
        if(!std::isfinite(coords[i]))
            Rcpp::stop("coords must contain only finite values");
    }
    if(flat_pair_edge_offsets.size() != pair_graph_distance.size() + 1)
        Rcpp::stop("flat_pair_edge_offsets must have length length(pair_graph_distance) + 1");
    if(flat_edge_u.size() != flat_edge_v.size() ||
       flat_edge_u.size() != flat_edge_coeff.size())
        Rcpp::stop("flat edge arrays must have the same length");
    if(flat_pair_edge_offsets[0] != 0)
        Rcpp::stop("flat_pair_edge_offsets must start at 0");
    if(flat_pair_edge_offsets[flat_pair_edge_offsets.size() - 1] != flat_edge_u.size())
        Rcpp::stop("flat_pair_edge_offsets must end at length(flat_edge_u)");
    for(int i = 0; i < flat_edge_u.size(); i++){
        if(flat_edge_u[i] < 0 || flat_edge_u[i] >= coords.nrow() ||
           flat_edge_v[i] < 0 || flat_edge_v[i] >= coords.nrow())
            Rcpp::stop("flat edge arrays must use 0-based vertex ids within [0, nrow(coords) - 1]");
        if(!std::isfinite(flat_edge_coeff[i]) || flat_edge_coeff[i] < 0.0)
            Rcpp::stop("flat_edge_coeff must contain finite values >= 0");
    }
    for(int i = 0; i < pair_graph_distance.size(); i++){
        if(!std::isfinite(pair_graph_distance[i]) || pair_graph_distance[i] <= 0.0)
            Rcpp::stop("pair_graph_distance must contain finite values > 0");
    }
    if(graph_edge_u.size() != graph_edge_v.size() ||
       graph_edge_u.size() != graph_edge_target.size())
        Rcpp::stop("graph edge arrays must have the same length");
    for(int i = 0; i < graph_edge_u.size(); i++){
        if(graph_edge_u[i] < 0 || graph_edge_u[i] >= coords.nrow() ||
           graph_edge_v[i] < 0 || graph_edge_v[i] >= coords.nrow())
            Rcpp::stop("graph edge arrays must use 0-based vertex ids within [0, nrow(coords) - 1]");
        if(!std::isfinite(graph_edge_target[i]) || graph_edge_target[i] <= 0.0)
            Rcpp::stop("graph_edge_target must contain finite values > 0");
    }
    if(repulsion_u.size() != repulsion_v.size() ||
       repulsion_u.size() != repulsion_target.size())
        Rcpp::stop("repulsion arrays must have the same length");
    for(int i = 0; i < repulsion_u.size(); i++){
        if(repulsion_u[i] < 0 || repulsion_u[i] >= coords.nrow() ||
           repulsion_v[i] < 0 || repulsion_v[i] >= coords.nrow())
            Rcpp::stop("repulsion arrays must use 0-based vertex ids within [0, nrow(coords) - 1]");
        if(!std::isfinite(repulsion_target[i]) || repulsion_target[i] < 0.0)
            Rcpp::stop("repulsion_target must contain finite values >= 0");
    }

    validate_geodesic_mds_args(max_iter,
                               edge_length_epsilon,
                               initial_step,
                               step_shrink,
                               armijo_factor,
                               grad_tol,
                               min_step);

    FlatGeodesicCacheView cache{
        Rcpp::as<std::vector<int>>(flat_pair_edge_offsets),
        Rcpp::as<std::vector<int>>(flat_edge_u),
        Rcpp::as<std::vector<int>>(flat_edge_v),
        Rcpp::as<std::vector<double>>(flat_edge_coeff),
        Rcpp::as<std::vector<double>>(pair_graph_distance)
    };
    std::vector<Point<>> current = matrix_to_points(coords);
    std::vector<Point<>> anchor;
    bool useAnchor = anchor_coords.isNotNull();
    if(useAnchor){
        Rcpp::NumericMatrix anchorMat = anchor_coords.get();
        if(anchorMat.nrow() != coords.nrow() || anchorMat.ncol() != coords.ncol())
            Rcpp::stop("anchor_coords must have the same dimensions as coords");
        for(int i = 0; i < anchorMat.size(); i++){
            if(!std::isfinite(anchorMat[i]))
                Rcpp::stop("anchor_coords must contain only finite values");
        }
        anchor = matrix_to_points(anchorMat);
    }
    std::vector<double> anchorSchedule = resolve_anchor_schedule(anchor_weights, max_iter);
    FlatSmoothnessView smoothness = build_flat_smoothness_view(
        smooth_adj_offsets,
        smooth_adj_vertices,
        coords.nrow()
    );
    std::vector<double> smoothSchedule = resolve_anchor_schedule(smooth_weights, max_iter);
    FlatEdgePenaltyView edgePenalty{
        Rcpp::as<std::vector<int>>(graph_edge_u),
        Rcpp::as<std::vector<int>>(graph_edge_v),
        Rcpp::as<std::vector<double>>(graph_edge_target)
    };
    std::vector<double> edgeSpringSchedule = resolve_anchor_schedule(edge_spring_weights, max_iter);
    FlatRepulsionView repulsion{
        Rcpp::as<std::vector<int>>(repulsion_u),
        Rcpp::as<std::vector<int>>(repulsion_v),
        Rcpp::as<std::vector<double>>(repulsion_target)
    };
    std::vector<double> repulsionSchedule = resolve_anchor_schedule(repulsion_weights, max_iter);
    int resolvedThreads = resolve_gmds_thread_count(n_threads, cache.pairCount());

    std::vector<int> trace_iteration;
    std::vector<double> trace_energy;
    std::vector<double> trace_gmds_energy;
    std::vector<double> trace_anchor_energy;
    std::vector<double> trace_edge_spring_energy;
    std::vector<double> trace_repulsion_energy;
    std::vector<int> trace_repulsion_pair_count;
    std::vector<int> trace_repulsion_active_pair_count;
    std::vector<double> trace_smooth_energy;
    std::vector<double> trace_gradient_norm;
    std::vector<double> trace_step;
    std::vector<bool> trace_accepted;
    std::vector<double> trace_anchor_weight;
    std::vector<double> trace_edge_spring_weight;
    std::vector<double> trace_repulsion_weight;
    std::vector<double> trace_smooth_weight;
    std::vector<std::vector<Point<>>> accepted_frames;
    accepted_frames.push_back(current);

    const double eps2 = edge_length_epsilon * edge_length_epsilon;
    const double gradTol2 = grad_tol * grad_tol;

    GeodesicMdsState state = evaluate_flat_state(
        current,
        cache,
        smoothness,
        edgePenalty,
        repulsion,
        eps2,
        useAnchor ? &anchor : nullptr,
        anchorSchedule[0],
        smoothSchedule[0],
        edgeSpringSchedule[0],
        repulsionSchedule[0],
        resolvedThreads
    );
    trace_iteration.push_back(0);
    trace_energy.push_back(state.energy);
    trace_gmds_energy.push_back(state.gmdsEnergy);
    trace_anchor_energy.push_back(state.anchorEnergy);
    trace_edge_spring_energy.push_back(state.edgeSpringEnergy);
    trace_repulsion_energy.push_back(state.repulsionEnergy);
    trace_repulsion_pair_count.push_back(state.repulsionPairCount);
    trace_repulsion_active_pair_count.push_back(state.repulsionActivePairCount);
    trace_smooth_energy.push_back(state.smoothEnergy);
    trace_gradient_norm.push_back(std::sqrt(state.gradNorm2));
    trace_step.push_back(NA_REAL);
    trace_accepted.push_back(true);
    trace_anchor_weight.push_back(anchorSchedule[0]);
    trace_edge_spring_weight.push_back(edgeSpringSchedule[0]);
    trace_repulsion_weight.push_back(repulsionSchedule[0]);
    trace_smooth_weight.push_back(smoothSchedule[0]);

    for(int iter = 1; iter <= max_iter; iter++){
        double iterAnchorWeight = anchorSchedule[static_cast<size_t>(iter)];
        double iterSmoothWeight = smoothSchedule[static_cast<size_t>(iter)];
        double iterEdgeSpringWeight = edgeSpringSchedule[static_cast<size_t>(iter)];
        double iterRepulsionWeight = repulsionSchedule[static_cast<size_t>(iter)];
        state = evaluate_flat_state(
            current,
            cache,
            smoothness,
            edgePenalty,
            repulsion,
            eps2,
            useAnchor ? &anchor : nullptr,
            iterAnchorWeight,
            iterSmoothWeight,
            iterEdgeSpringWeight,
            iterRepulsionWeight,
            resolvedThreads
        );
        if(!std::isfinite(state.energy) || state.gradNorm2 <= gradTol2)
            break;

        double step = initial_step;
        bool accepted = false;
        std::vector<Point<>> proposal = current;
        GeodesicMdsState candidate = state;

        while(std::isfinite(step) && step >= min_step){
            proposal = current;
            for(size_t i = 0; i < proposal.size(); i++)
                proposal[i] -= state.gradient[i] * step;
            if(recenter)
                recenter_points(proposal);

            candidate = evaluate_flat_state(
                proposal,
                cache,
                smoothness,
                edgePenalty,
                repulsion,
                eps2,
                useAnchor ? &anchor : nullptr,
                iterAnchorWeight,
                iterSmoothWeight,
                iterEdgeSpringWeight,
                iterRepulsionWeight,
                resolvedThreads
            );
            double targetEnergy = state.energy - armijo_factor * step * state.gradNorm2;
            if(std::isfinite(candidate.energy) && candidate.energy <= targetEnergy){
                accepted = true;
                break;
            }
            step *= step_shrink;
        }

        trace_iteration.push_back(iter);
        trace_energy.push_back(accepted ? candidate.energy : state.energy);
        trace_gmds_energy.push_back(accepted ? candidate.gmdsEnergy : state.gmdsEnergy);
        trace_anchor_energy.push_back(accepted ? candidate.anchorEnergy : state.anchorEnergy);
        trace_edge_spring_energy.push_back(accepted ? candidate.edgeSpringEnergy : state.edgeSpringEnergy);
        trace_repulsion_energy.push_back(accepted ? candidate.repulsionEnergy : state.repulsionEnergy);
        trace_repulsion_pair_count.push_back(accepted ? candidate.repulsionPairCount : state.repulsionPairCount);
        trace_repulsion_active_pair_count.push_back(accepted ? candidate.repulsionActivePairCount : state.repulsionActivePairCount);
        trace_smooth_energy.push_back(accepted ? candidate.smoothEnergy : state.smoothEnergy);
        trace_gradient_norm.push_back(std::sqrt(accepted ? candidate.gradNorm2 : state.gradNorm2));
        trace_step.push_back(accepted ? step : NA_REAL);
        trace_accepted.push_back(accepted);
        trace_anchor_weight.push_back(iterAnchorWeight);
        trace_edge_spring_weight.push_back(iterEdgeSpringWeight);
        trace_repulsion_weight.push_back(iterRepulsionWeight);
        trace_smooth_weight.push_back(iterSmoothWeight);

        if(!accepted)
            break;

        current = proposal;
        state = candidate;
        if(return_trace)
            accepted_frames.push_back(current);
    }

    if(!return_trace){
        accepted_frames.clear();
        accepted_frames.push_back(current);
    }

    return Rcpp::List::create(
        Rcpp::_["coords"] = points_to_matrix(current, coords.ncol()),
        Rcpp::_["trace"] = build_trace_df(trace_iteration,
                                          trace_energy,
                                          trace_gmds_energy,
                                          trace_anchor_energy,
                                          trace_edge_spring_energy,
                                          trace_repulsion_energy,
                                          trace_repulsion_pair_count,
                                          trace_repulsion_active_pair_count,
                                          trace_smooth_energy,
                                          trace_gradient_norm,
                                          trace_step,
                                          trace_accepted,
                                          trace_anchor_weight,
                                          trace_edge_spring_weight,
                                          trace_repulsion_weight,
                                          trace_smooth_weight),
        Rcpp::_["frames"] = build_frame_list(accepted_frames, coords.ncol()),
        Rcpp::_["final_anchor_weight"] = trace_anchor_weight.back(),
        Rcpp::_["final_edge_spring_weight"] = trace_edge_spring_weight.back(),
        Rcpp::_["final_repulsion_weight"] = trace_repulsion_weight.back(),
        Rcpp::_["final_smoothness_weight"] = trace_smooth_weight.back(),
        Rcpp::_["n_threads_used"] = resolvedThreads
    );
}

// [[Rcpp::export]]
Rcpp::List grip_geodesic_mds_flat_repulsion_stats_cpp(
    Rcpp::IntegerVector repulsion_u,
    Rcpp::IntegerVector repulsion_v,
    Rcpp::NumericVector repulsion_target,
    Rcpp::NumericMatrix coords,
    double edge_length_epsilon = 1e-8,
    double repulsion_weight = 1.0)
{
    if(coords.ncol() != 2 && coords.ncol() != 3)
        Rcpp::stop("coords must have 2 or 3 columns");
    if(coords.nrow() <= 0)
        Rcpp::stop("coords must have at least one row");
    for(int i = 0; i < coords.size(); i++){
        if(!std::isfinite(coords[i]))
            Rcpp::stop("coords must contain only finite values");
    }
    if(!std::isfinite(edge_length_epsilon) || edge_length_epsilon < 0.0)
        Rcpp::stop("edge_length_epsilon must be finite and >= 0");
    if(!std::isfinite(repulsion_weight) || repulsion_weight < 0.0)
        Rcpp::stop("repulsion_weight must be finite and >= 0");
    if(repulsion_u.size() != repulsion_v.size() ||
       repulsion_u.size() != repulsion_target.size())
        Rcpp::stop("repulsion arrays must have the same length");
    for(int i = 0; i < repulsion_u.size(); i++){
        if(repulsion_u[i] < 0 || repulsion_u[i] >= coords.nrow() ||
           repulsion_v[i] < 0 || repulsion_v[i] >= coords.nrow())
            Rcpp::stop("repulsion arrays must use 0-based vertex ids within [0, nrow(coords) - 1]");
        if(!std::isfinite(repulsion_target[i]) || repulsion_target[i] < 0.0)
            Rcpp::stop("repulsion_target must contain finite values >= 0");
    }

    FlatRepulsionView repulsion{
        Rcpp::as<std::vector<int>>(repulsion_u),
        Rcpp::as<std::vector<int>>(repulsion_v),
        Rcpp::as<std::vector<double>>(repulsion_target)
    };
    std::vector<Point<>> current = matrix_to_points(coords);
    std::vector<Point<>> gradient(current.size());
    for(size_t i = 0; i < gradient.size(); i++)
        gradient[i].set_to_zero();
    double energy = 0.0;
    int activePairCount = 0;
    accumulate_flat_repulsion(current,
                              repulsion,
                              edge_length_epsilon * edge_length_epsilon,
                              repulsion_weight,
                              gradient,
                              energy,
                              activePairCount);

    return Rcpp::List::create(
        Rcpp::_["repulsion_weight"] = repulsion_weight,
        Rcpp::_["energy"] = energy,
        Rcpp::_["gradient"] = points_to_matrix(gradient, coords.ncol()),
        Rcpp::_["pair_count"] = static_cast<int>(repulsion.pairCount()),
        Rcpp::_["active_pair_count"] = activePairCount
    );
}
