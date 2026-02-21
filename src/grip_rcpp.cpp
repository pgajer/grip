#include <Rcpp.h>
#include <ctime>
#include <cmath>

#include "Graph.h"
#include "DrawGraph.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix grip_layout_cpp(Rcpp::IntegerMatrix edges,
                                    Rcpp::Nullable<Rcpp::NumericVector> edge_weights,
                                    int n,
                                    int dim,
                                    std::string engine,
                                    std::string placement,
                                    int rounds,
                                    int final_rounds,
                                    int num_init,
                                    int num_nbrs,
                                    double r,
                                    double s,
                                    int tinit_factor,
                                    Rcpp::Nullable<int> seed)
{
    if(dim != 2 && dim != 3)
        Rcpp::stop("dim must be 2 or 3");
    if(edges.ncol() != 2)
        Rcpp::stop("edges must have exactly 2 columns");
    if(n <= 0)
        Rcpp::stop("n must be positive");

    if(num_init <= 0)
        num_init = 1;
    if(num_nbrs <= 0)
        num_nbrs = 1;
    if(rounds <= 0)
        rounds = 1;
    if(final_rounds <= 0)
        final_rounds = 1;
    if(tinit_factor <= 0)
        tinit_factor = 1;

    std::vector<std::pair<size_tt, size_tt>> edge_list;
    edge_list.reserve(edges.nrow());
    for(int i = 0; i < edges.nrow(); i++){
        int u = edges(i, 0);
        int v = edges(i, 1);
        if(u <= 0 || v <= 0)
            Rcpp::stop("edges must be 1-based vertex ids");
        edge_list.emplace_back(static_cast<size_tt>(u - 1),
                               static_cast<size_tt>(v - 1));
    }

    std::vector<coord_t> weights_vec;
    std::vector<coord_t> *weights_ptr = nullptr;
    if(edge_weights.isNotNull()){
        Rcpp::NumericVector w = edge_weights.get();
        if(w.size() != edges.nrow())
            Rcpp::stop("edge_weights length must match number of edges");
        weights_vec.reserve(w.size());
        for(int i = 0; i < w.size(); i++){
            double wi = w[i];
            if(!std::isfinite(wi) || wi <= 0.0)
                Rcpp::stop("edge_weights must contain finite values > 0; invalid value %.17g at edge index %d",
                           wi,
                           i + 1);
            weights_vec.push_back(static_cast<coord_t>(wi));
        }
        weights_ptr = &weights_vec;
    }

    Graph graph;
    unsigned int seed_val = seed.isNotNull()
        ? static_cast<unsigned int>(Rcpp::as<int>(seed))
        : static_cast<unsigned int>(std::time(nullptr));
    graph.sfast_Rand(seed_val);
    graph.from_edge_list(static_cast<size_tt>(n), edge_list, weights_ptr);

    size_tt engf = (engine == "mish_v6") ? 13 : 12;
    size_tt placement_mode =
        (placement == "circle") ? PLACEMENT_CIRCLE : PLACEMENT_BARYCENTER;

    DrawGraph dg(graph,
                 static_cast<size_tt>(dim),
                 static_cast<size_tt>(rounds),
                 static_cast<size_tt>(final_rounds),
                 static_cast<size_tt>(tinit_factor),
                 engf,
                 static_cast<size_tt>(num_init),
                 static_cast<size_tt>(num_nbrs),
                 r,
                 s,
                 placement_mode,
                 false);

    if(engf == 12)
        dg.mish_engine_v5();
    else
        dg.mish_engine_v6();

    Rcpp::NumericMatrix out(n, dim);
    Point<> *pos = dg.get_Pos();
    for(int i = 0; i < n; i++){
        out(i, 0) = pos[i].getX();
        if(dim > 1) out(i, 1) = pos[i].getY();
        if(dim > 2) out(i, 2) = pos[i].getZ();
    }

    return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix grip_layout_adj_cpp(Rcpp::List adj_list,
                                        Rcpp::Nullable<Rcpp::List> weight_list,
                                        int n,
                                        int dim,
                                        std::string engine,
                                        std::string placement,
                                        int rounds,
                                        int final_rounds,
                                        int num_init,
                                        int num_nbrs,
                                        double r,
                                        double s,
                                        int tinit_factor,
                                        Rcpp::Nullable<int> seed)
{
    if(dim != 2 && dim != 3)
        Rcpp::stop("dim must be 2 or 3");
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

    if(num_init <= 0)
        num_init = 1;
    if(num_nbrs <= 0)
        num_nbrs = 1;
    if(rounds <= 0)
        rounds = 1;
    if(final_rounds <= 0)
        final_rounds = 1;
    if(tinit_factor <= 0)
        tinit_factor = 1;

    std::vector<std::vector<size_tt>> adj(n);
    std::vector<std::vector<coord_t>> weights;
    if(useWeights)
        weights.resize(n);
    for(int i = 0; i < n; i++){
        Rcpp::IntegerVector neigh = adj_list[i];
        adj[i].reserve(neigh.size());
        Rcpp::NumericVector w;
        if(useWeights){
            w = weight_list_val[i];
            if(neigh.size() != w.size())
                Rcpp::stop("weight_list must be parallel to adj_list");
            weights[i].reserve(w.size());
        }
        for(int j = 0; j < neigh.size(); j++){
            int v = neigh[j];
            if(v <= 0 || v > n)
                Rcpp::stop("adj_list must be 1-based and within [1, n]");
            adj[i].push_back(static_cast<size_tt>(v - 1));
            if(useWeights){
                double wj = w[j];
                if(!std::isfinite(wj) || wj <= 0.0)
                    Rcpp::stop("weight_list must contain finite values > 0; invalid value %.17g at weight_list[[%d]][%d]",
                               wj,
                               i + 1,
                               j + 1);
                weights[i].push_back(static_cast<coord_t>(wj));
            }
        }
    }

    Graph graph;
    unsigned int seed_val = seed.isNotNull()
        ? static_cast<unsigned int>(Rcpp::as<int>(seed))
        : static_cast<unsigned int>(std::time(nullptr));
    graph.sfast_Rand(seed_val);
    graph.from_adj_list(static_cast<size_tt>(n), adj, useWeights ? &weights : nullptr);

    size_tt engf = (engine == "mish_v6") ? 13 : 12;
    size_tt placement_mode =
        (placement == "circle") ? PLACEMENT_CIRCLE : PLACEMENT_BARYCENTER;

    DrawGraph dg(graph,
                 static_cast<size_tt>(dim),
                 static_cast<size_tt>(rounds),
                 static_cast<size_tt>(final_rounds),
                 static_cast<size_tt>(tinit_factor),
                 engf,
                 static_cast<size_tt>(num_init),
                 static_cast<size_tt>(num_nbrs),
                 r,
                 s,
                 placement_mode,
                 false);

    if(engf == 12)
        dg.mish_engine_v5();
    else
        dg.mish_engine_v6();

    Rcpp::NumericMatrix out(n, dim);
    Point<> *pos = dg.get_Pos();
    for(int i = 0; i < n; i++){
        out(i, 0) = pos[i].getX();
        if(dim > 1) out(i, 1) = pos[i].getY();
        if(dim > 2) out(i, 2) = pos[i].getZ();
    }

    return out;
}
