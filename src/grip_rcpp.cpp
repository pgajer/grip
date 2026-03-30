#include <Rcpp.h>
#include <algorithm>
#include <ctime>
#include <cmath>

#include "Graph.h"
#include "DrawGraph.h"

namespace {

void validate_tuning_args(int num_nbrs,
                          double r,
                          double s,
                          double repulsion_factor)
{
    if(num_nbrs <= 0)
        Rcpp::stop("num_nbrs must be a positive integer");
    if(!std::isfinite(r) || r < 0.0 || r > 1.0)
        Rcpp::stop("r must be finite and in [0, 1]");
    if(!std::isfinite(s) || s < 0.0)
        Rcpp::stop("s must be finite and >= 0");
    if(!std::isfinite(repulsion_factor) || repulsion_factor < 0.0)
        Rcpp::stop("repulsion_factor must be finite and >= 0");
}

void validate_globalrep_tuning_args(int num_nbrs,
                                    double r,
                                    double s,
                                    double repulsion_factor,
                                    double coarse_repulsion_factor,
                                    int coarse_repulsion_sample,
                                    int coarse_repulsion_exact_below,
                                    double final_anchor_factor,
                                    double final_move_scale_after_first)
{
    validate_tuning_args(num_nbrs, r, s, repulsion_factor);
    if(!std::isfinite(coarse_repulsion_factor) || coarse_repulsion_factor < 0.0)
        Rcpp::stop("coarse_repulsion_factor must be finite and >= 0");
    if(coarse_repulsion_sample <= 0)
        Rcpp::stop("coarse_repulsion_sample must be a positive integer");
    if(coarse_repulsion_exact_below <= 0)
        Rcpp::stop("coarse_repulsion_exact_below must be a positive integer");
    if(!std::isfinite(final_anchor_factor) || final_anchor_factor < 0.0)
        Rcpp::stop("final_anchor_factor must be finite and >= 0");
    if(!std::isfinite(final_move_scale_after_first) ||
       final_move_scale_after_first < 0.0 ||
       final_move_scale_after_first > 1.0)
        Rcpp::stop("final_move_scale_after_first must be finite and in [0, 1]");
}

size_tt final_stage_mode_from_string(const std::string &final_mode)
{
    if(final_mode == "fr")
        return FINAL_STAGE_FR;
    if(final_mode == "kk_repulse")
        return FINAL_STAGE_KK_REPULSE;
    Rcpp::stop("final_mode must be either 'fr' or 'kk_repulse'");
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix grip_layout_cpp(Rcpp::IntegerMatrix edges,
                                    Rcpp::Nullable<Rcpp::NumericVector> edge_weights,
                                    int n,
                                    int dim,
                                    std::string placement,
                                    int rounds,
                                    int final_rounds,
                                    int num_init,
                                    int num_nbrs,
                                    double r,
                                    double s,
                                    double repulsion_factor,
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
    validate_tuning_args(num_nbrs, r, s, repulsion_factor);
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

    size_tt placement_mode =
        (placement == "circle") ? PLACEMENT_CIRCLE : PLACEMENT_BARYCENTER;

    DrawGraph dg(graph,
                 static_cast<size_tt>(dim),
                 static_cast<size_tt>(rounds),
                 static_cast<size_tt>(final_rounds),
                 static_cast<size_tt>(tinit_factor),
                 static_cast<size_tt>(num_init),
                 static_cast<size_tt>(num_nbrs),
                 r,
                 s,
                 repulsion_factor,
                 placement_mode,
                 false);

    dg.mish_engine();

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
                                        std::string placement,
                                        int rounds,
                                        int final_rounds,
                                        int num_init,
                                        int num_nbrs,
                                        double r,
                                        double s,
                                        double repulsion_factor,
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
    validate_tuning_args(num_nbrs, r, s, repulsion_factor);
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

    size_tt placement_mode =
        (placement == "circle") ? PLACEMENT_CIRCLE : PLACEMENT_BARYCENTER;

    DrawGraph dg(graph,
                 static_cast<size_tt>(dim),
                 static_cast<size_tt>(rounds),
                 static_cast<size_tt>(final_rounds),
                 static_cast<size_tt>(tinit_factor),
                 static_cast<size_tt>(num_init),
                 static_cast<size_tt>(num_nbrs),
                 r,
                 s,
                 repulsion_factor,
                 placement_mode,
                 false);

    dg.mish_engine();

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
Rcpp::NumericMatrix grip_layout_globalrep_adj_cpp(
    Rcpp::List adj_list,
    Rcpp::Nullable<Rcpp::List> weight_list,
    int n,
    int dim,
    std::string placement,
    int rounds,
    int final_rounds,
    int num_init,
    int num_nbrs,
    double r,
    double s,
    double repulsion_factor,
    double coarse_repulsion_factor,
    int coarse_repulsion_sample,
    int coarse_repulsion_exact_below,
    double final_anchor_factor,
    double final_move_scale_after_first,
    std::string final_mode,
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
    validate_globalrep_tuning_args(num_nbrs,
                                   r,
                                   s,
                                   repulsion_factor,
                                   coarse_repulsion_factor,
                                   coarse_repulsion_sample,
                                   coarse_repulsion_exact_below,
                                   final_anchor_factor,
                                   final_move_scale_after_first);
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

    size_tt placement_mode =
        (placement == "circle") ? PLACEMENT_CIRCLE : PLACEMENT_BARYCENTER;
    size_tt final_stage_mode = final_stage_mode_from_string(final_mode);

    DrawGraph dg(graph,
                 static_cast<size_tt>(dim),
                 static_cast<size_tt>(rounds),
                 static_cast<size_tt>(final_rounds),
                 static_cast<size_tt>(tinit_factor),
                 static_cast<size_tt>(num_init),
                 static_cast<size_tt>(num_nbrs),
                 r,
                 s,
                 repulsion_factor,
                 placement_mode,
                 false,
                 final_stage_mode,
                 coarse_repulsion_factor,
                 static_cast<size_tt>(coarse_repulsion_sample),
                 static_cast<size_tt>(coarse_repulsion_exact_below),
                 final_anchor_factor,
                 final_move_scale_after_first);

    dg.mish_engine();

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
Rcpp::List grip_layout_trace_adj_cpp(Rcpp::List adj_list,
                                     Rcpp::Nullable<Rcpp::List> weight_list,
                                     int n,
                                     int dim,
                                     std::string placement,
                                     int rounds,
                                     int final_rounds,
                                     int num_init,
                                     int num_nbrs,
                                     double r,
                                     double s,
                                     double repulsion_factor,
                                     int tinit_factor,
                                     Rcpp::Nullable<int> seed,
                                     std::string trace,
                                     int trace_every)
{
    if(dim != 2 && dim != 3)
        Rcpp::stop("dim must be 2 or 3");
    if(n <= 0)
        Rcpp::stop("n must be positive");
    if(adj_list.size() != n)
        Rcpp::stop("adj_list length must match n");
    if(trace != "round" && trace != "level")
        Rcpp::stop("trace must be either 'round' or 'level'");
    if(trace_every <= 0)
        Rcpp::stop("trace_every must be a positive integer");

    bool useWeights = weight_list.isNotNull();
    Rcpp::List weight_list_val;
    if(useWeights){
        weight_list_val = weight_list.get();
        if(weight_list_val.size() != n)
            Rcpp::stop("weight_list length must match n");
    }

    if(num_init <= 0)
        num_init = 1;
    validate_tuning_args(num_nbrs, r, s, repulsion_factor);
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

    size_tt placement_mode =
        (placement == "circle") ? PLACEMENT_CIRCLE : PLACEMENT_BARYCENTER;

    DrawGraph dg(graph,
                 static_cast<size_tt>(dim),
                 static_cast<size_tt>(rounds),
                 static_cast<size_tt>(final_rounds),
                 static_cast<size_tt>(tinit_factor),
                 static_cast<size_tt>(num_init),
                 static_cast<size_tt>(num_nbrs),
                 r,
                 s,
                 repulsion_factor,
                 placement_mode,
                 false);
    dg.configure_trace(trace == "round" ? TRACE_ROUND : TRACE_LEVEL,
                       static_cast<size_tt>(trace_every));

    dg.mish_engine();

    Rcpp::NumericMatrix out(n, dim);
    Point<> *pos = dg.get_Pos();
    for(int i = 0; i < n; i++){
        out(i, 0) = pos[i].getX();
        if(dim > 1) out(i, 1) = pos[i].getY();
        if(dim > 2) out(i, 2) = pos[i].getZ();
    }

    const auto &trace_frames = dg.get_trace_frames();
    Rcpp::List frames(trace_frames.size());
    for(size_t i = 0; i < trace_frames.size(); i++){
        Rcpp::NumericMatrix frame(n, dim);
        std::copy(trace_frames[i].begin(), trace_frames[i].end(), frame.begin());
        frames[i] = frame;
    }

    Rcpp::IntegerVector frame_ids(trace_frames.size());
    for(int i = 0; i < frame_ids.size(); i++)
        frame_ids[i] = i + 1;

    Rcpp::DataFrame meta = Rcpp::DataFrame::create(
        Rcpp::_["frame"] = frame_ids,
        Rcpp::_["phase"] = Rcpp::wrap(dg.get_trace_phases()),
        Rcpp::_["level_index"] = Rcpp::wrap(dg.get_trace_level_indices()),
        Rcpp::_["misf_level"] = Rcpp::wrap(dg.get_trace_misf_levels()),
        Rcpp::_["round_in_level"] = Rcpp::wrap(dg.get_trace_rounds()),
        Rcpp::_["active_vertices"] = Rcpp::wrap(dg.get_trace_active_counts()),
        Rcpp::_["stringsAsFactors"] = false
    );

    return Rcpp::List::create(
        Rcpp::_["final"] = out,
        Rcpp::_["frames"] = frames,
        Rcpp::_["meta"] = meta
    );
}

// [[Rcpp::export]]
Rcpp::List grip_layout_globalrep_trace_adj_cpp(Rcpp::List adj_list,
                                               Rcpp::Nullable<Rcpp::List> weight_list,
                                               int n,
                                               int dim,
                                               std::string placement,
                                               int rounds,
                                               int final_rounds,
                                               int num_init,
                                               int num_nbrs,
                                               double r,
                                               double s,
                                               double repulsion_factor,
                                               double coarse_repulsion_factor,
                                               int coarse_repulsion_sample,
                                               int coarse_repulsion_exact_below,
                                               double final_anchor_factor,
                                               double final_move_scale_after_first,
                                               std::string final_mode,
                                               int tinit_factor,
                                               Rcpp::Nullable<int> seed,
                                               std::string trace,
                                               int trace_every)
{
    if(dim != 2 && dim != 3)
        Rcpp::stop("dim must be 2 or 3");
    if(n <= 0)
        Rcpp::stop("n must be positive");
    if(adj_list.size() != n)
        Rcpp::stop("adj_list length must match n");
    if(trace != "round" && trace != "level")
        Rcpp::stop("trace must be either 'round' or 'level'");
    if(trace_every <= 0)
        Rcpp::stop("trace_every must be a positive integer");

    bool useWeights = weight_list.isNotNull();
    Rcpp::List weight_list_val;
    if(useWeights){
        weight_list_val = weight_list.get();
        if(weight_list_val.size() != n)
            Rcpp::stop("weight_list length must match n");
    }

    if(num_init <= 0)
        num_init = 1;
    validate_globalrep_tuning_args(num_nbrs,
                                   r,
                                   s,
                                   repulsion_factor,
                                   coarse_repulsion_factor,
                                   coarse_repulsion_sample,
                                   coarse_repulsion_exact_below,
                                   final_anchor_factor,
                                   final_move_scale_after_first);
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

    size_tt placement_mode =
        (placement == "circle") ? PLACEMENT_CIRCLE : PLACEMENT_BARYCENTER;
    size_tt final_stage_mode = final_stage_mode_from_string(final_mode);

    DrawGraph dg(graph,
                 static_cast<size_tt>(dim),
                 static_cast<size_tt>(rounds),
                 static_cast<size_tt>(final_rounds),
                 static_cast<size_tt>(tinit_factor),
                 static_cast<size_tt>(num_init),
                 static_cast<size_tt>(num_nbrs),
                 r,
                 s,
                 repulsion_factor,
                 placement_mode,
                 false,
                 final_stage_mode,
                 coarse_repulsion_factor,
                 static_cast<size_tt>(coarse_repulsion_sample),
                 static_cast<size_tt>(coarse_repulsion_exact_below),
                 final_anchor_factor,
                 final_move_scale_after_first);
    dg.configure_trace(trace == "round" ? TRACE_ROUND : TRACE_LEVEL,
                       static_cast<size_tt>(trace_every));

    dg.mish_engine();

    Rcpp::NumericMatrix out(n, dim);
    Point<> *pos = dg.get_Pos();
    for(int i = 0; i < n; i++){
        out(i, 0) = pos[i].getX();
        if(dim > 1) out(i, 1) = pos[i].getY();
        if(dim > 2) out(i, 2) = pos[i].getZ();
    }

    const auto &trace_frames = dg.get_trace_frames();
    Rcpp::List frames(trace_frames.size());
    for(size_t i = 0; i < trace_frames.size(); i++){
        Rcpp::NumericMatrix frame(n, dim);
        std::copy(trace_frames[i].begin(), trace_frames[i].end(), frame.begin());
        frames[i] = frame;
    }

    Rcpp::IntegerVector frame_ids(trace_frames.size());
    for(int i = 0; i < frame_ids.size(); i++)
        frame_ids[i] = i + 1;

    Rcpp::DataFrame meta = Rcpp::DataFrame::create(
        Rcpp::_["frame"] = frame_ids,
        Rcpp::_["phase"] = Rcpp::wrap(dg.get_trace_phases()),
        Rcpp::_["level_index"] = Rcpp::wrap(dg.get_trace_level_indices()),
        Rcpp::_["misf_level"] = Rcpp::wrap(dg.get_trace_misf_levels()),
        Rcpp::_["round_in_level"] = Rcpp::wrap(dg.get_trace_rounds()),
        Rcpp::_["active_vertices"] = Rcpp::wrap(dg.get_trace_active_counts()),
        Rcpp::_["stringsAsFactors"] = false
    );

    return Rcpp::List::create(
        Rcpp::_["final"] = out,
        Rcpp::_["frames"] = frames,
        Rcpp::_["meta"] = meta
    );
}
