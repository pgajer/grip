#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

#include "Graph.h"
#include "Point.h"

namespace {

struct GeodesicPathEdge {
    size_tt u;
    size_tt v;
};

struct GeodesicPairCache {
    coord_t graphDistance;
    std::vector<GeodesicPathEdge> pathEdges;
};

struct GeodesicMdsState {
    double energy;
    double gradNorm2;
    std::vector<Point<>> gradient;
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
                    GeodesicPathEdge{pathVertices[edgeIndex - 1], pathVertices[edgeIndex]}
                );
            }
            pairs.push_back(std::move(pair));
        }
    }

    return pairs;
}

GeodesicMdsState evaluate_state(const std::vector<Point<>> &coords,
                                const std::vector<GeodesicPairCache> &pairs,
                                double eps2)
{
    GeodesicMdsState state;
    state.energy = 0.0;
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
            h += len;
        }

        double resid = h - pair.graphDistance;
        state.energy += 0.5 * resid * resid;
        for(size_t edgeIndex = 0; edgeIndex < pair.pathEdges.size(); edgeIndex++){
            if(edgeLens[edgeIndex] <= 0.0)
                continue;
            const GeodesicPathEdge &edgeRef = pair.pathEdges[edgeIndex];
            Point<> stepVec = edgeDiffs[edgeIndex] * (resid / edgeLens[edgeIndex]);
            state.gradient[static_cast<size_t>(edgeRef.u)] += stepVec;
            state.gradient[static_cast<size_t>(edgeRef.v)] -= stepVec;
        }
    }

    for(size_t i = 0; i < state.gradient.size(); i++)
        state.gradNorm2 += state.gradient[i].fnorm2();

    return state;
}

Rcpp::DataFrame build_trace_df(const std::vector<int> &iteration,
                               const std::vector<double> &energy,
                               const std::vector<double> &gradient_norm,
                               const std::vector<double> &step,
                               const std::vector<bool> &accepted)
{
    return Rcpp::DataFrame::create(
        Rcpp::_["iteration"] = iteration,
        Rcpp::_["energy"] = energy,
        Rcpp::_["gradient_norm"] = gradient_norm,
        Rcpp::_["step"] = step,
        Rcpp::_["accepted"] = accepted,
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
    std::vector<double> trace_gradient_norm;
    std::vector<double> trace_step;
    std::vector<bool> trace_accepted;
    std::vector<std::vector<Point<>>> accepted_frames;
    accepted_frames.push_back(current);

    const double eps2 = edge_length_epsilon * edge_length_epsilon;
    const double gradTol2 = grad_tol * grad_tol;

    GeodesicMdsState state = evaluate_state(current, pairs, eps2);
    trace_iteration.push_back(0);
    trace_energy.push_back(state.energy);
    trace_gradient_norm.push_back(std::sqrt(state.gradNorm2));
    trace_step.push_back(NA_REAL);
    trace_accepted.push_back(true);

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

            candidate = evaluate_state(proposal, pairs, eps2);
            double targetEnergy = state.energy - armijo_factor * step * state.gradNorm2;
            if(std::isfinite(candidate.energy) && candidate.energy <= targetEnergy){
                accepted = true;
                break;
            }
            step *= step_shrink;
        }

        trace_iteration.push_back(iter);
        trace_energy.push_back(accepted ? candidate.energy : state.energy);
        trace_gradient_norm.push_back(std::sqrt(accepted ? candidate.gradNorm2 : state.gradNorm2));
        trace_step.push_back(accepted ? step : NA_REAL);
        trace_accepted.push_back(accepted);

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
                                          trace_gradient_norm,
                                          trace_step,
                                          trace_accepted),
        Rcpp::_["frames"] = build_frame_list(accepted_frames, coords.ncol())
    );
}
