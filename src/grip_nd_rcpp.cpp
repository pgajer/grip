#include <Rcpp.h>

#include <cmath>
#include <ctime>
#include <limits>
#include <queue>

#include "DrawGraphND.h"
#include "GraphND.h"
#include "WeightedMisfND.h"

namespace {

void validate_weighted_nd_tuning(int dim,
                                 int rounds,
                                 int final_rounds,
                                 int num_init,
                                 int num_nbrs,
                                 double r,
                                 double s,
                                 double repulsion_factor,
                                 int tinit_factor)
{
    if(dim < 2)
        Rcpp::stop("dim must be >= 2");
    if(rounds < 0)
        Rcpp::stop("rounds must be a non-negative integer");
    if(final_rounds < 0)
        Rcpp::stop("final_rounds must be a non-negative integer");
    if(num_init < dim + 1)
        Rcpp::stop("num_init must be at least dim + 1 for ND weighted layout");
    if(num_nbrs <= 0)
        Rcpp::stop("num_nbrs must be a positive integer");
    if(!std::isfinite(r) || r < 0.0 || r > 1.0)
        Rcpp::stop("r must be finite and in [0, 1]");
    if(!std::isfinite(s) || s < 0.0)
        Rcpp::stop("s must be finite and >= 0");
    if(!std::isfinite(repulsion_factor) || repulsion_factor < 0.0)
        Rcpp::stop("repulsion_factor must be finite and >= 0");
    if(tinit_factor <= 0)
        Rcpp::stop("tinit_factor must be a positive integer");
}

int parse_placement_mode_nd(const std::string &placement, int dim)
{
    if(placement == "barycenter")
        return gripnd::PLACEMENT_BARYCENTER_ND;
    if(placement == "circle"){
        if(dim != 2)
            Rcpp::stop("placement = 'circle' requires dim = 2");
        return gripnd::PLACEMENT_CIRCLE_ND;
    }
    Rcpp::stop("placement must be 'barycenter' or 'circle'");
}

int parse_insertion_anchor_scope_nd(const std::string &scope)
{
    if(scope == "any_higher")
        return gripnd::INSERT_ANCHOR_SCOPE_ANY_HIGHER_ND;
    if(scope == "prev_misf")
        return gripnd::INSERT_ANCHOR_SCOPE_PREV_MISF_ND;
    Rcpp::stop("insertion_anchor_scope must be 'any_higher' or 'prev_misf'");
}

int parse_insertion_anchor_strategy_nd(const std::string &strategy)
{
    if(strategy == "first")
        return gripnd::INSERT_ANCHOR_STRATEGY_FIRST_ND;
    if(strategy == "distance_band")
        return gripnd::INSERT_ANCHOR_STRATEGY_DISTANCE_BAND_ND;
    if(strategy == "balanced_band")
        return gripnd::INSERT_ANCHOR_STRATEGY_BALANCED_BAND_ND;
    if(strategy == "spread_prev")
        return gripnd::INSERT_ANCHOR_STRATEGY_SPREAD_PREV_ND;
    Rcpp::stop("insertion_anchor_strategy must be 'first', 'distance_band', 'balanced_band', or 'spread_prev'");
}

int parse_level0_insertion_mode_nd(const std::string &mode)
{
    if(mode == "inherit")
        return gripnd::LEVEL0_INSERT_INHERIT_ND;
    if(mode == "barycenter")
        return gripnd::LEVEL0_INSERT_BARYCENTER_ND;
    if(mode == "least_squares")
        return gripnd::LEVEL0_INSERT_LEAST_SQUARES_ND;
    Rcpp::stop("level0_insertion_mode must be 'inherit', 'barycenter', or 'least_squares'");
}

int parse_final_stage_mode_nd(const std::string &mode)
{
    if(mode == "fr")
        return gripnd::FINAL_STAGE_FR_ND;
    if(mode == "kk_repulse")
        return gripnd::FINAL_STAGE_KK_REPULSE_ND;
    Rcpp::stop("final_mode must be 'fr' or 'kk_repulse'");
}

gripnd::GraphND graph_from_r_lists(Rcpp::List adj_list,
                                   Rcpp::List weight_list,
                                   int n)
{
    if(adj_list.size() != n)
        Rcpp::stop("adj_list length must equal n");
    if(weight_list.size() != n)
        Rcpp::stop("weight_list length must equal n");

    gripnd::GraphND graph(static_cast<std::size_t>(n));
    for(int i = 0; i < n; i++){
        Rcpp::IntegerVector nbr = adj_list[i];
        Rcpp::NumericVector wt = weight_list[i];
        if(nbr.size() != wt.size())
            Rcpp::stop("weight_list entries must be parallel to adj_list entries");

        for(int k = 0; k < nbr.size(); k++){
            const int v = nbr[k];
            const double w = wt[k];
            if(v < 1 || v > n)
                Rcpp::stop("adj_list contains a vertex outside [1, n]");
            if(v == i + 1)
                continue;
            if(!std::isfinite(w) || w <= 0.0)
                Rcpp::stop("weight_list must contain finite values > 0");
            graph.add_directed_edge(
                static_cast<gripnd::vertex_t>(i),
                static_cast<gripnd::vertex_t>(v - 1),
                w
            );
        }
    }
    return graph;
}

Rcpp::NumericMatrix points_to_matrix(const std::vector<gripnd::PointND> &coords,
                                     int n,
                                     int dim)
{
    Rcpp::NumericMatrix out(n, dim);
    for(int i = 0; i < n; i++){
        for(int d = 0; d < dim; d++)
            out(i, d) = coords[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
    }
    return out;
}

Rcpp::List trace_to_list(const gripnd::LayoutTraceND &trace,
                         int n,
                         int dim)
{
    Rcpp::List frames(trace.frames.size());
    for(std::size_t f = 0; f < trace.frames.size(); f++){
        Rcpp::NumericMatrix frame(n, dim);
        const std::vector<double> &flat = trace.frames[f];
        for(int i = 0; i < n; i++){
            for(int d = 0; d < dim; d++)
                frame(i, d) = flat[static_cast<std::size_t>(i * dim + d)];
        }
        frames[f] = frame;
    }

    return Rcpp::List::create(
        Rcpp::Named("frames") = frames,
        Rcpp::Named("phase") = trace.phases,
        Rcpp::Named("round") = trace.rounds,
        Rcpp::Named("level_index") = trace.level_indices,
        Rcpp::Named("misf_level") = trace.misf_levels,
        Rcpp::Named("active_vertices") = trace.active_counts
    );
}

Rcpp::NumericMatrix nested_vectors_to_matrix(
    const std::vector<std::vector<double>> &values,
    int dim)
{
    Rcpp::NumericMatrix out(values.size(), dim);
    for(std::size_t i = 0; i < values.size(); i++){
        for(int d = 0; d < dim; d++)
            out(i, d) = values[i][static_cast<std::size_t>(d)];
    }
    return out;
}

Rcpp::List refinement_step_trace_to_list(
    const gripnd::RefinementStepTraceND &trace,
    int dim)
{
    return Rcpp::List::create(
        Rcpp::Named("trace") = Rcpp::DataFrame::create(
            Rcpp::Named("level_index") = trace.level_indices,
            Rcpp::Named("misf_level") = trace.misf_levels,
            Rcpp::Named("round_in_level") = trace.rounds,
            Rcpp::Named("active_vertices") = trace.active_counts,
            Rcpp::Named("order_index") = trace.order_indices,
            Rcpp::Named("vertex") = trace.vertices,
            Rcpp::Named("heat_before") = trace.heat_before,
            Rcpp::Named("heat_after") = trace.heat_after,
            Rcpp::Named("old_cos_before") = trace.old_cos_before,
            Rcpp::Named("old_cos_after") = trace.old_cos_after,
            Rcpp::Named("old_disp_norm_before") = trace.old_disp_norm_before,
            Rcpp::Named("pre_temp_disp_norm") = trace.pre_temp_disp_norm,
            Rcpp::Named("attraction_edges") = trace.attraction_edges,
            Rcpp::Named("repulsion_neighbors") = trace.repulsion_neighbors,
            Rcpp::Named("stringsAsFactors") = false
        ),
        Rcpp::Named("coords_before") = nested_vectors_to_matrix(trace.coords_before, dim),
        Rcpp::Named("pre_temp_disp") = nested_vectors_to_matrix(trace.pre_temp_disp, dim),
        Rcpp::Named("attraction_disp") = nested_vectors_to_matrix(trace.attraction_disp, dim),
        Rcpp::Named("repulsion_disp") = nested_vectors_to_matrix(trace.repulsion_disp, dim),
        Rcpp::Named("applied_disp") = nested_vectors_to_matrix(trace.applied_disp, dim),
        Rcpp::Named("coords_after") = nested_vectors_to_matrix(trace.coords_after, dim),
        Rcpp::Named("attraction_terms") = Rcpp::DataFrame::create(
            Rcpp::Named("parent_row") = trace.attraction_term_parent_rows,
            Rcpp::Named("term_index") = trace.attraction_term_indices,
            Rcpp::Named("vertex") = trace.attraction_term_vertices,
            Rcpp::Named("neighbor") = trace.attraction_term_neighbors,
            Rcpp::Named("weight") = trace.attraction_term_weights,
            Rcpp::Named("norm2") = trace.attraction_term_norm2,
            Rcpp::Named("desired") = trace.attraction_term_desired,
            Rcpp::Named("desired2") = trace.attraction_term_desired2,
            Rcpp::Named("scale") = trace.attraction_term_scale,
            Rcpp::Named("stringsAsFactors") = false
        ),
        Rcpp::Named("attraction_term_delta") =
            nested_vectors_to_matrix(trace.attraction_term_delta, dim),
        Rcpp::Named("attraction_term_step") =
            nested_vectors_to_matrix(trace.attraction_term_step, dim),
        Rcpp::Named("attraction_term_cumulative") =
            nested_vectors_to_matrix(trace.attraction_term_cumulative, dim)
    );
}

Rcpp::List insertion_trace_to_list(const gripnd::InsertionTraceND &trace,
                                   int dim)
{
    return Rcpp::List::create(
        Rcpp::Named("trace") = Rcpp::DataFrame::create(
            Rcpp::Named("level_index") = trace.level_indices,
            Rcpp::Named("misf_level") = trace.misf_levels,
            Rcpp::Named("previous_active_vertices") = trace.previous_active_counts,
            Rcpp::Named("active_vertices") = trace.active_counts,
            Rcpp::Named("order_index") = trace.order_indices,
            Rcpp::Named("vertex") = trace.vertices,
            Rcpp::Named("root_depth") = trace.root_depths,
            Rcpp::Named("anchor_count_requested") = trace.anchor_count_requested,
            Rcpp::Named("anchor_count_used") = trace.anchor_count_used,
            Rcpp::Named("insertion_mode") = trace.insertion_modes,
            Rcpp::Named("local_kk_steps") = trace.local_kk_steps,
            Rcpp::Named("anchors") = trace.anchors,
            Rcpp::Named("old_disp_norm_initial") = trace.old_disp_norm_initial,
            Rcpp::Named("old_disp_norm_after") = trace.old_disp_norm_after,
            Rcpp::Named("stringsAsFactors") = false
        ),
        Rcpp::Named("coords_initial") =
            nested_vectors_to_matrix(trace.coords_initial, dim),
        Rcpp::Named("coords_after") =
            nested_vectors_to_matrix(trace.coords_after, dim),
        Rcpp::Named("old_disp_initial") =
            nested_vectors_to_matrix(trace.old_disp_initial, dim),
        Rcpp::Named("old_disp_after") =
            nested_vectors_to_matrix(trace.old_disp_after, dim)
    );
}

int strict_misf_prefix_weighted_nd(const gripnd::GraphND &graph,
                                   const std::vector<gripnd::vertex_t> &order,
                                   int top_size,
                                   int level)
{
    if(level == 0 || top_size < 2)
        return top_size;

    const double radius = (level <= 1)
        ? 1.0
        : std::pow(2.0, static_cast<double>(level - 1));
    const double tol = 1e-10;
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<int> in_prefix(graph.size(), -1);

    for(int idx = 0; idx < top_size; idx++){
        const gripnd::vertex_t root = order[static_cast<std::size_t>(idx)];
        if(idx > 0){
            struct Node {
                double dist;
                gripnd::vertex_t vert;
            };
            struct Greater {
                bool operator()(const Node &lhs, const Node &rhs) const
                {
                    if(lhs.dist != rhs.dist)
                        return lhs.dist > rhs.dist;
                    return lhs.vert > rhs.vert;
                }
            };

            std::vector<double> dist(graph.size(), inf);
            std::priority_queue<Node, std::vector<Node>, Greater> pq;
            dist[root] = 0.0;
            pq.push(Node{0.0, root});

            while(!pq.empty()){
                Node node = pq.top();
                pq.pop();
                if(node.dist > dist[node.vert] + tol)
                    continue;
                if(node.dist > radius + tol)
                    break;

                if(node.vert != root &&
                   in_prefix[node.vert] >= 0 &&
                   in_prefix[node.vert] < idx)
                    return idx;

                for(const gripnd::NeighborND &nbr : graph.neighbors(node.vert)){
                    double alt = node.dist + nbr.weight;
                    if(alt > radius + tol)
                        continue;
                    if(alt + tol < dist[nbr.vertex]){
                        dist[nbr.vertex] = alt;
                        pq.push(Node{alt, nbr.vertex});
                    }
                }
            }
        }
        in_prefix[root] = idx;
    }

    return top_size;
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix grip_layout_weighted_nd_adj_cpp(Rcpp::List adj_list,
                                                    Rcpp::List weight_list,
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
                                                    double final_anchor_factor,
                                                    double final_move_scale_after_first,
                                                    std::string final_mode,
                                                    int metric_neighbor_cap,
                                                    int insertion_anchor_count,
                                                    std::string insertion_anchor_scope,
                                                    std::string insertion_anchor_strategy,
                                                    std::string level0_insertion_mode,
                                                    int level0_anchor_count,
                                                    int level0_local_kk_steps,
                                                    Rcpp::Nullable<int> seed = R_NilValue)
{
    if(n <= 0)
        Rcpp::stop("n must be a positive integer");
    validate_weighted_nd_tuning(
        dim,
        rounds,
        final_rounds,
        num_init,
        num_nbrs,
        r,
        s,
        repulsion_factor,
        tinit_factor
    );
    if(insertion_anchor_count <= 0)
        Rcpp::stop("insertion_anchor_count must be a positive integer");
    if(level0_anchor_count <= 0)
        Rcpp::stop("level0_anchor_count must be a positive integer");
    if(level0_local_kk_steps < 0)
        Rcpp::stop("level0_local_kk_steps must be a non-negative integer");
    if(!std::isfinite(final_move_scale_after_first) ||
       final_move_scale_after_first < 0.0 ||
       final_move_scale_after_first > 1.0)
        Rcpp::stop("final_move_scale_after_first must be in [0, 1]");
    if(!std::isfinite(final_anchor_factor) || final_anchor_factor < 0.0)
        Rcpp::stop("final_anchor_factor must be finite and >= 0");
    if(metric_neighbor_cap < 0)
        Rcpp::stop("metric_neighbor_cap must be >= 0");

    const unsigned int seed_value = seed.isNotNull()
        ? static_cast<unsigned int>(Rcpp::as<int>(seed))
        : static_cast<unsigned int>(std::time(NULL));
    const int placement_mode = parse_placement_mode_nd(placement, dim);
    const int insertion_anchor_scope_code =
        parse_insertion_anchor_scope_nd(insertion_anchor_scope);
    const int insertion_anchor_strategy_code =
        parse_insertion_anchor_strategy_nd(insertion_anchor_strategy);
    const int level0_insertion_mode_code =
        parse_level0_insertion_mode_nd(level0_insertion_mode);
    const int final_stage_mode_code = parse_final_stage_mode_nd(final_mode);

    gripnd::GraphND graph = graph_from_r_lists(adj_list, weight_list, n);
    gripnd::DrawGraphND drawer(
        graph,
        dim,
        rounds,
        final_rounds,
        num_init,
        num_nbrs,
        r,
        s,
        repulsion_factor,
        tinit_factor,
        final_anchor_factor,
        final_move_scale_after_first,
        final_stage_mode_code,
        metric_neighbor_cap,
        placement_mode,
        insertion_anchor_count,
        insertion_anchor_scope_code,
        insertion_anchor_strategy_code,
        level0_insertion_mode_code,
        level0_anchor_count,
        level0_local_kk_steps,
        seed_value
    );

    return points_to_matrix(drawer.layout(), n, dim);
}

// [[Rcpp::export]]
Rcpp::List grip_layout_weighted_nd_trace_adj_cpp(Rcpp::List adj_list,
                                                 Rcpp::List weight_list,
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
                                                 double final_anchor_factor,
                                                 double final_move_scale_after_first,
                                                 std::string final_mode,
                                                 int metric_neighbor_cap,
                                                 int insertion_anchor_count,
                                                 std::string insertion_anchor_scope,
                                                 std::string insertion_anchor_strategy,
                                                 std::string level0_insertion_mode,
                                                 int level0_anchor_count,
                                                 int level0_local_kk_steps,
                                                 Rcpp::Nullable<int> seed = R_NilValue,
                                                 int trace_every = 1,
                                                 bool refinement_step_trace = false,
                                                 int refinement_step_level_index = -1,
                                                 int refinement_step_misf_level = -1,
                                                 int refinement_step_round_start = -1,
                                                 int refinement_step_round_end = -1)
{
    if(n <= 0)
        Rcpp::stop("n must be a positive integer");
    validate_weighted_nd_tuning(
        dim,
        rounds,
        final_rounds,
        num_init,
        num_nbrs,
        r,
        s,
        repulsion_factor,
        tinit_factor
    );
    if(trace_every <= 0)
        Rcpp::stop("trace_every must be a positive integer");
    if(insertion_anchor_count <= 0)
        Rcpp::stop("insertion_anchor_count must be a positive integer");
    if(level0_anchor_count <= 0)
        Rcpp::stop("level0_anchor_count must be a positive integer");
    if(level0_local_kk_steps < 0)
        Rcpp::stop("level0_local_kk_steps must be a non-negative integer");
    if(!std::isfinite(final_move_scale_after_first) ||
       final_move_scale_after_first < 0.0 ||
       final_move_scale_after_first > 1.0)
        Rcpp::stop("final_move_scale_after_first must be in [0, 1]");
    if(!std::isfinite(final_anchor_factor) || final_anchor_factor < 0.0)
        Rcpp::stop("final_anchor_factor must be finite and >= 0");
    if(metric_neighbor_cap < 0)
        Rcpp::stop("metric_neighbor_cap must be >= 0");

    const unsigned int seed_value = seed.isNotNull()
        ? static_cast<unsigned int>(Rcpp::as<int>(seed))
        : static_cast<unsigned int>(std::time(NULL));
    const int placement_mode = parse_placement_mode_nd(placement, dim);
    const int insertion_anchor_scope_code =
        parse_insertion_anchor_scope_nd(insertion_anchor_scope);
    const int insertion_anchor_strategy_code =
        parse_insertion_anchor_strategy_nd(insertion_anchor_strategy);
    const int level0_insertion_mode_code =
        parse_level0_insertion_mode_nd(level0_insertion_mode);
    const int final_stage_mode_code = parse_final_stage_mode_nd(final_mode);

    gripnd::GraphND graph = graph_from_r_lists(adj_list, weight_list, n);
    gripnd::DrawGraphND drawer(
        graph,
        dim,
        rounds,
        final_rounds,
        num_init,
        num_nbrs,
        r,
        s,
        repulsion_factor,
        tinit_factor,
        final_anchor_factor,
        final_move_scale_after_first,
        final_stage_mode_code,
        metric_neighbor_cap,
        placement_mode,
        insertion_anchor_count,
        insertion_anchor_scope_code,
        insertion_anchor_strategy_code,
        level0_insertion_mode_code,
        level0_anchor_count,
        level0_local_kk_steps,
        seed_value
    );

    gripnd::LayoutTraceND trace;
    if(refinement_step_trace){
        drawer.configure_refinement_step_trace(
            refinement_step_level_index,
            refinement_step_misf_level,
            refinement_step_round_start,
            refinement_step_round_end
        );
        drawer.configure_insertion_trace(true);
    }
    std::vector<gripnd::PointND> coords = drawer.layout(&trace, trace_every);
    Rcpp::List trace_list = trace_to_list(trace, n, dim);
    Rcpp::List step_trace = refinement_step_trace_to_list(
        drawer.get_refinement_step_trace(),
        dim
    );
    return Rcpp::List::create(
        Rcpp::Named("final") = points_to_matrix(coords, n, dim),
        Rcpp::Named("frames") = trace_list["frames"],
        Rcpp::Named("phase") = trace_list["phase"],
        Rcpp::Named("round") = trace_list["round"],
        Rcpp::Named("level_index") = trace_list["level_index"],
        Rcpp::Named("misf_level") = trace_list["misf_level"],
        Rcpp::Named("active_vertices") = trace_list["active_vertices"],
        Rcpp::Named("refinement_step_trace") = step_trace,
        Rcpp::Named("insertion_trace") =
            insertion_trace_to_list(drawer.get_insertion_trace(), dim)
    );
}

// [[Rcpp::export]]
Rcpp::List grip_build_weighted_misf_nd_adj_cpp(Rcpp::List adj_list,
                                               Rcpp::List weight_list,
                                               int n,
                                               int num_init,
                                               int num_nbrs,
                                               Rcpp::Nullable<int> seed = R_NilValue)
{
    if(n <= 0)
        Rcpp::stop("n must be a positive integer");
    if(num_init <= 0)
        Rcpp::stop("num_init must be a positive integer");
    if(num_nbrs <= 0)
        Rcpp::stop("num_nbrs must be a positive integer");

    const unsigned long seed_value = seed.isNotNull()
        ? static_cast<unsigned long>(Rcpp::as<int>(seed))
        : static_cast<unsigned long>(std::time(NULL));

    gripnd::GraphND graph = graph_from_r_lists(adj_list, weight_list, n);
    gripnd::WeightedMisfND misf = gripnd::build_weighted_misf_nd(
        graph,
        num_init,
        num_nbrs,
        seed_value
    );

    std::vector<int> level_size = misf.level_size;
    std::vector<int> vertex_depth_values = misf.vertex_depth;
    if(level_size.size() >= 2){
        const int top_level = static_cast<int>(level_size.size()) - 1;
        const int strict_size = strict_misf_prefix_weighted_nd(
            graph,
            misf.order,
            level_size[static_cast<std::size_t>(top_level)],
            top_level
        );
        if(strict_size < level_size[static_cast<std::size_t>(top_level)]){
            for(int idx = strict_size;
                idx < level_size[static_cast<std::size_t>(top_level)];
                idx++){
                vertex_depth_values[misf.order[static_cast<std::size_t>(idx)]] = top_level - 1;
            }
            level_size[static_cast<std::size_t>(top_level)] = strict_size;
        }
    }

    Rcpp::List levels(level_size.size());
    Rcpp::CharacterVector level_names(level_size.size());
    Rcpp::IntegerVector misf_size(level_size.size());
    Rcpp::IntegerVector num_nbrs_schedule(misf.num_nbrs_schedule.size());
    for(std::size_t level = 0; level < level_size.size(); level++){
        const int size = level_size[level];
        Rcpp::IntegerVector verts(size);
        for(int i = 0; i < size; i++)
            verts[i] = static_cast<int>(misf.order[static_cast<std::size_t>(i)]) + 1;
        levels[level] = verts;
        level_names[level] = "V" + std::to_string(level);
        misf_size[level] = size;
        num_nbrs_schedule[level] = misf.num_nbrs_schedule[level];
    }
    levels.attr("names") = level_names;

    Rcpp::IntegerVector vertex_depth_r(n);
    Rcpp::IntegerVector mish_order(n);
    for(int i = 0; i < n; i++){
        vertex_depth_r[i] = vertex_depth_values[static_cast<std::size_t>(i)];
        mish_order[i] = static_cast<int>(misf.order[static_cast<std::size_t>(i)]) + 1;
    }

    return Rcpp::List::create(
        Rcpp::Named("levels") = levels,
        Rcpp::Named("vertex_depth") = vertex_depth_r,
        Rcpp::Named("mish_order") = mish_order,
        Rcpp::Named("misf_size") = misf_size,
        Rcpp::Named("num_nbrs_schedule") = num_nbrs_schedule,
        Rcpp::Named("misf_height") = misf.height,
        Rcpp::Named("top_level_size") = misf_size[misf_size.size() - 1]
    );
}
