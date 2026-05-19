#include "DrawGraphND.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <limits>
#include <queue>
#include <sstream>
#include <unordered_set>

namespace gripnd {
namespace {

constexpr double kLegacyEdgeND = 32.0;
constexpr double kMetricTolND = 1e-10;

double dot_point_nd(const PointND &lhs, const PointND &rhs)
{
    double out = 0.0;
    for(std::size_t d = 0; d < lhs.dim(); d++)
        out += lhs[d] * rhs[d];
    return out;
}

double round_legacy_norm_nd(double value)
{
    return value > 0.0 ? std::floor(value + 0.5) : -std::floor(0.5 - value);
}

} // namespace

DrawGraphND::DrawGraphND(const GraphND &graph,
                         int dim,
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
                         int tinit_factor,
                         double final_anchor_factor,
                         double final_move_scale_after_first,
                         int final_stage_mode,
                         int metric_neighbor_cap,
                         int placement_mode,
                         int insertion_anchor_count,
                         int insertion_anchor_scope,
                         int insertion_anchor_strategy,
                         int level0_insertion_mode,
                         int level0_anchor_count,
                         int level0_local_kk_steps,
                         int lgkk_multiscale_rounds,
                         int lgkk_rounds_coarse,
                         int lgkk_rounds_pre_final,
                         int lgkk_rounds_final,
                         int lgkk_local_nbrs,
                         int lgkk_landmark_count,
                         int lgkk_scope,
                         int lgkk_active_limit,
                         unsigned int seed)
    : graph_(graph),
      dim_(dim),
      rounds_(rounds),
      final_rounds_(final_rounds),
      num_init_(num_init),
      num_nbrs_(num_nbrs),
      r_(r),
      s_(s),
      repulsion_factor_(repulsion_factor),
      coarse_repulsion_factor_(std::max(0.0, coarse_repulsion_factor)),
      coarse_repulsion_sample_(std::max(0, coarse_repulsion_sample)),
      coarse_repulsion_exact_below_(std::max(0, coarse_repulsion_exact_below)),
      tinit_factor_(tinit_factor),
      final_anchor_factor_(std::max(0.0, final_anchor_factor)),
      final_move_scale_after_first_(
          std::max(0.0, std::min(1.0, final_move_scale_after_first))),
      final_stage_mode_(final_stage_mode == FINAL_STAGE_KK_REPULSE_ND
                            ? FINAL_STAGE_KK_REPULSE_ND
                            : FINAL_STAGE_FR_ND),
      metric_neighbor_cap_(std::max(0, metric_neighbor_cap)),
      placement_mode_(placement_mode == PLACEMENT_CIRCLE_ND
                          ? PLACEMENT_CIRCLE_ND
                          : PLACEMENT_BARYCENTER_ND),
      insertion_anchor_count_(std::max(1, insertion_anchor_count)),
      insertion_anchor_scope_(insertion_anchor_scope == INSERT_ANCHOR_SCOPE_PREV_MISF_ND
                                  ? INSERT_ANCHOR_SCOPE_PREV_MISF_ND
                                  : INSERT_ANCHOR_SCOPE_ANY_HIGHER_ND),
      insertion_anchor_strategy_(
          insertion_anchor_strategy == INSERT_ANCHOR_STRATEGY_DISTANCE_BAND_ND
              ? INSERT_ANCHOR_STRATEGY_DISTANCE_BAND_ND
              : insertion_anchor_strategy == INSERT_ANCHOR_STRATEGY_BALANCED_BAND_ND
                    ? INSERT_ANCHOR_STRATEGY_BALANCED_BAND_ND
                    : insertion_anchor_strategy == INSERT_ANCHOR_STRATEGY_SPREAD_PREV_ND
                          ? INSERT_ANCHOR_STRATEGY_SPREAD_PREV_ND
                          : INSERT_ANCHOR_STRATEGY_FIRST_ND),
      level0_insertion_mode_(
          level0_insertion_mode == LEVEL0_INSERT_BARYCENTER_ND
              ? LEVEL0_INSERT_BARYCENTER_ND
              : level0_insertion_mode == LEVEL0_INSERT_LEAST_SQUARES_ND
                    ? LEVEL0_INSERT_LEAST_SQUARES_ND
                    : LEVEL0_INSERT_INHERIT_ND),
      level0_anchor_count_(std::max(1, level0_anchor_count)),
      level0_local_kk_steps_(std::max(0, level0_local_kk_steps)),
      lgkk_multiscale_rounds_(std::max(0, lgkk_multiscale_rounds)),
      lgkk_rounds_coarse_(std::max(0, lgkk_rounds_coarse)),
      lgkk_rounds_pre_final_(std::max(0, lgkk_rounds_pre_final)),
      lgkk_rounds_final_(std::max(0, lgkk_rounds_final)),
      lgkk_local_nbrs_(std::max(0, lgkk_local_nbrs)),
      lgkk_landmark_count_(std::max(0, lgkk_landmark_count)),
      lgkk_scope_(lgkk_scope == LGKK_SCOPE_COARSE_ND
                      ? LGKK_SCOPE_COARSE_ND
                      : LGKK_SCOPE_ALL_ND),
      lgkk_active_limit_(std::max(1, lgkk_active_limit)),
      seed_(seed),
      rng_(seed),
      legacy_rng_state_(seed),
      box_size_(1.0),
      box2_size_(3),
      coords_(graph.size(), PointND(static_cast<std::size_t>(dim))),
      disp_(graph.size(), PointND(static_cast<std::size_t>(dim))),
      old_disp_(graph.size(), PointND(static_cast<std::size_t>(dim))),
      disp_norm_(graph.size(), 0.0),
      old_disp_norm_(graph.size(), 0.0),
      heat_(graph.size(), kLegacyEdgeND / std::max(1, tinit_factor)),
      old_cos_(graph.size(), 1.0),
      final_anchor_ready_(false),
      final_anchor_pos_(graph.size(), PointND(static_cast<std::size_t>(dim))),
      metric_neighbors_cache_(graph.size()),
      metric_neighbors_cached_(graph.size(), 0),
      lgkk_cache_active_count_(0),
      lgkk_cache_misf_level_(-1),
      lgkk_cache_scale_l0_(1.0),
      refinement_step_trace_enabled_(false),
      refinement_step_trace_level_index_(-1),
      refinement_step_trace_misf_level_(-1),
      refinement_step_trace_round_start_(-1),
      refinement_step_trace_round_end_(-1),
      insertion_trace_enabled_(false),
      last_attraction_disp_(static_cast<std::size_t>(dim), 0.0),
      last_repulsion_disp_(static_cast<std::size_t>(dim), 0.0)
{
}

std::vector<PointND> DrawGraphND::layout()
{
    return layout(nullptr, 1);
}

std::vector<PointND> DrawGraphND::layout(LayoutTraceND *trace, int trace_every)
{
    trace_every = std::max(1, trace_every);
    const WeightedMisfND misf = build_weighted_misf_nd(
        graph_,
        num_init_,
        num_nbrs_,
        seed_
    );
    initialize_multiscale_trace(misf);

    int level_index = 1;
    int misf_level = misf.height;
    int active_count = misf.num_init;
    int previous_active_count = active_count;
    int current_level_rounds = rounds_;
    record_trace(trace, "init", 0, active_count, level_index, misf_level);
    if(misf_level == 0 && active_count == static_cast<int>(graph_.size()))
        prepare_final_anchors(active_count);
    refine_legacy_weighted_level(misf,
                                 current_level_rounds,
                                 active_count,
                                 trace,
                                 trace_every,
                                 level_index,
                                 misf_level);
    current_level_rounds = lgkk_refine_level(misf,
                                             active_count,
                                             misf_level,
                                             current_level_rounds,
                                             trace,
                                             trace_every,
                                             level_index);

    while(active_count < static_cast<int>(graph_.size())){
        level_index++;
        if(misf_level > 0 &&
           misf.level_size[static_cast<std::size_t>(misf_level)] !=
               static_cast<int>(graph_.size())){
            misf_level--;
            active_count = misf.level_size[static_cast<std::size_t>(misf_level)];
        } else {
            active_count = static_cast<int>(graph_.size());
        }
        reset_active_heat(previous_active_count);
        insert_level_vertices(misf,
                              previous_active_count,
                              active_count,
                              level_index,
                              misf_level);
        current_level_rounds = scheduled_rounds(active_count);
        record_trace(trace, "level_start", 0, active_count, level_index, misf_level);
        if(misf_level == 0 && active_count == static_cast<int>(graph_.size()))
            prepare_final_anchors(active_count);
        refine_legacy_weighted_level(misf,
                                     current_level_rounds,
                                     active_count,
                                     trace,
                                     trace_every,
                                     level_index,
                                     misf_level);
        current_level_rounds = lgkk_refine_level(misf,
                                                 active_count,
                                                 misf_level,
                                                 current_level_rounds,
                                                 trace,
                                                 trace_every,
                                                 level_index);
        previous_active_count = active_count;
    }
    record_trace(trace,
                 "final",
                 current_level_rounds,
                 active_count,
                 level_index,
                 misf_level);
    return coords_;
}

void DrawGraphND::configure_refinement_step_trace(int level_index,
                                                  int misf_level,
                                                  int round_start,
                                                  int round_end)
{
    refinement_step_trace_enabled_ = true;
    refinement_step_trace_level_index_ = level_index;
    refinement_step_trace_misf_level_ = misf_level;
    refinement_step_trace_round_start_ = round_start;
    refinement_step_trace_round_end_ = std::max(round_start, round_end);
    refinement_step_trace_ = RefinementStepTraceND();
}

void DrawGraphND::configure_insertion_trace(bool enabled)
{
    insertion_trace_enabled_ = enabled;
    insertion_trace_ = InsertionTraceND();
}

void DrawGraphND::initialize()
{
    const double base_scale = graph_.median_edge_weight();
    const double spread =
        base_scale * std::max(1.0, std::sqrt(static_cast<double>(num_init_))) *
        std::max(1.0, static_cast<double>(tinit_factor_) / 2.0);

    std::uniform_real_distribution<double> unif(-spread, spread);
    for(PointND &p : coords_){
        for(int d = 0; d < dim_; d++)
            p[static_cast<std::size_t>(d)] = unif(rng_);
    }
}

void DrawGraphND::refine(const std::vector<EdgeND> &edges, int total_rounds)
{
    if(total_rounds <= 0 || coords_.empty())
        return;

    const double target_scale = graph_.median_edge_weight();
    const double base_move = std::max(1e-4, target_scale * std::max(0.01, r_));
    const double edge_rate = 0.35 / std::max(1.0, std::sqrt(static_cast<double>(dim_)));
    const double repulse_base =
        repulsion_factor_ * 0.006 * std::max(1.0, s_) /
        std::max(1.0, std::sqrt(static_cast<double>(dim_)));

    for(int round = 0; round < total_rounds; round++){
        const double progress = static_cast<double>(round) /
            static_cast<double>(std::max(1, total_rounds - 1));
        const double cooling = 1.0 - 0.85 * progress;
        std::vector<PointND> disp(coords_.size(), PointND(static_cast<std::size_t>(dim_)));

        add_edge_forces(edges, disp, edge_rate);
        if(repulsion_factor_ > 0.0)
            add_repulsion_forces(disp, repulse_base * cooling);

        apply_displacements(disp, base_move * cooling);
    }
}

void DrawGraphND::refine_trace_level(const std::vector<EdgeND> &edges,
                                     int level_rounds,
                                     int active_count,
                                     LayoutTraceND *trace,
                                     int trace_every,
                                     int level_index,
                                     int misf_level)
{
    if(level_rounds <= 0 || coords_.empty())
        return;

    const double target_scale = graph_.median_edge_weight();
    const double base_move = std::max(1e-4, target_scale * std::max(0.01, r_));
    const double edge_rate = 0.35 / std::max(1.0, std::sqrt(static_cast<double>(dim_)));
    const double repulse_base =
        repulsion_factor_ * 0.006 * std::max(1.0, s_) /
        std::max(1.0, std::sqrt(static_cast<double>(dim_)));
    std::vector<char> active(coords_.size(), 0);
    for(int i = 0; i < active_count && i < static_cast<int>(trace_order_.size()); i++)
        active[trace_order_[static_cast<std::size_t>(i)]] = 1;

    for(int round = 0; round < level_rounds; round++){
        const double progress = static_cast<double>(round) /
            static_cast<double>(std::max(1, level_rounds - 1));
        const double cooling = 1.0 - 0.85 * progress;
        std::vector<PointND> disp(coords_.size(), PointND(static_cast<std::size_t>(dim_)));
        std::vector<EdgeND> active_edges;
        active_edges.reserve(edges.size());
        for(const EdgeND &edge : edges){
            if(active[edge.source] && active[edge.target])
                active_edges.push_back(edge);
        }

        add_edge_forces(active_edges, disp, edge_rate);
        for(vertex_t i = 0; i < disp.size(); i++){
            if(!active[i])
                disp[i].fill(0.0);
        }
        if(repulsion_factor_ > 0.0 && active_count == static_cast<int>(coords_.size()))
            add_repulsion_forces(disp, repulse_base * cooling);

        apply_displacements(disp, base_move * cooling, active);
        if(trace != nullptr && ((round + 1) % trace_every) == 0)
            record_trace(trace, "round", round + 1, active_count, level_index, misf_level);
    }
}

void DrawGraphND::add_edge_forces(const std::vector<EdgeND> &edges,
                                  std::vector<PointND> &disp,
                                  double edge_rate) const
{
    const double eps = 1e-9;
    for(const EdgeND &edge : edges){
        const PointND &a = coords_[edge.source];
        const PointND &b = coords_[edge.target];
        const double dist = std::sqrt(std::max(eps, squared_distance(a, b)));
        const double target = std::max(eps, edge.weight);
        const double pull = edge_rate * (dist - target) / dist;

        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            const double delta = b[dd] - a[dd];
            const double step = pull * delta;
            disp[edge.source][dd] += step;
            disp[edge.target][dd] -= step;
        }
    }
}

void DrawGraphND::add_repulsion_forces(std::vector<PointND> &disp, double repulse_rate) const
{
    const double eps = 1e-6;
    const double cap = std::max(16, num_nbrs_) * graph_.median_edge_weight();
    const double cap2 = cap * cap;
    for(vertex_t i = 0; i < coords_.size(); i++){
        for(vertex_t j = i + 1; j < coords_.size(); j++){
            const double dist2 = std::max(eps, squared_distance(coords_[i], coords_[j]));
            if(dist2 > cap2)
                continue;
            const double rate = repulse_rate / dist2;
            for(int d = 0; d < dim_; d++){
                const std::size_t dd = static_cast<std::size_t>(d);
                double delta = coords_[i][dd] - coords_[j][dd];
                if(std::fabs(delta) < eps)
                    delta = (d == 0) ? eps : 0.0;
                const double step = rate * delta;
                disp[i][dd] += step;
                disp[j][dd] -= step;
            }
        }
    }
}

void DrawGraphND::apply_displacements(const std::vector<PointND> &disp, double max_move)
{
    for(vertex_t i = 0; i < coords_.size(); i++){
        const double len = disp[i].norm();
        const double scale = (len > max_move && len > 0.0) ? (max_move / len) : 1.0;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            coords_[i][dd] += disp[i][dd] * scale;
            if(!std::isfinite(coords_[i][dd]))
                coords_[i][dd] = 0.0;
        }
    }
}

void DrawGraphND::apply_displacements(const std::vector<PointND> &disp,
                                      double max_move,
                                      const std::vector<char> &active)
{
    for(vertex_t i = 0; i < coords_.size(); i++){
        if(i >= active.size() || !active[i])
            continue;
        const double len = disp[i].norm();
        const double scale = (len > max_move && len > 0.0) ? (max_move / len) : 1.0;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            coords_[i][dd] += disp[i][dd] * scale;
            if(!std::isfinite(coords_[i][dd]))
                coords_[i][dd] = 0.0;
        }
    }
}

void DrawGraphND::record_trace(LayoutTraceND *trace,
                               const std::string &phase,
                               int round,
                               int active_count,
                               int level_index,
                               int misf_level) const
{
    if(trace == nullptr)
        return;

    active_count = std::max(0, std::min(active_count, static_cast<int>(coords_.size())));
    std::vector<double> frame(coords_.size() * static_cast<std::size_t>(dim_),
                              std::numeric_limits<double>::quiet_NaN());
    for(int i = 0; i < active_count; i++){
        const vertex_t vert = trace_order_.empty()
            ? static_cast<vertex_t>(i)
            : trace_order_[static_cast<std::size_t>(i)];
        const PointND &point = coords_[vert];
        for(int d = 0; d < dim_; d++)
            frame[static_cast<std::size_t>(vert * dim_ + d)] =
                point[static_cast<std::size_t>(d)];
    }
    trace->frames.push_back(frame);
    trace->phases.push_back(phase);
    trace->rounds.push_back(round);
    trace->level_indices.push_back(level_index);
    trace->misf_levels.push_back(misf_level);
    trace->active_counts.push_back(active_count);
}

void DrawGraphND::reset_active_heat(int active_count)
{
    const double tinit = kLegacyEdgeND / std::max(1, tinit_factor_);
    for(int i = 0; i < active_count && i < static_cast<int>(trace_order_.size()); i++)
        heat_[trace_order_[static_cast<std::size_t>(i)]] = tinit;
}

void DrawGraphND::refine_legacy_weighted_level(const WeightedMisfND &misf,
                                               int level_rounds,
                                               int active_count,
                                               LayoutTraceND *trace,
                                               int trace_every,
                                               int level_index,
                                               int misf_level)
{
    if(level_rounds <= 0 || active_count <= 0)
        return;

    active_count = std::min(active_count, static_cast<int>(misf.order.size()));
    for(int round = 0; round < level_rounds; round++){
        const int round_in_level = round + 1;
        const bool record_steps =
            should_record_refinement_step(level_index, misf_level, round_in_level);
        std::vector<std::size_t> step_rows(static_cast<std::size_t>(active_count),
                                           static_cast<std::size_t>(-1));
        for(int i = 0; i < active_count; i++){
            const vertex_t vert = misf.order[static_cast<std::size_t>(i)];
            const PointND coord_before = coords_[vert];
            const double heat_before = heat_[vert];
            const double old_cos_before = old_cos_[vert];
            const double old_disp_norm_before = old_disp_norm_[vert];
            const std::vector<MetricNeighborND> neighbors =
                metric_neighbors_for_level(vert, misf, misf_level);
            if(misf_level == 0 &&
               final_stage_mode_ == FINAL_STAGE_KK_REPULSE_ND)
                legacy_weighted_kk_final_displacement(vert, neighbors, active_count);
            else if(misf_level == 0)
                legacy_weighted_fr_displacement(vert, neighbors, active_count);
            else
                legacy_weighted_kk_displacement(vert, neighbors, active_count, misf_level);

            const PointND pre_temp_disp = disp_[vert];
            const double pre_temp_disp_norm = disp_norm_[vert];
            update_local_temperature(vert);
            const double heat_after = heat_[vert];
            const double old_cos_after = old_cos_[vert];
            old_disp_[vert] = disp_[vert];
            old_disp_norm_[vert] = disp_norm_[vert];
            for(int d = 0; d < dim_; d++)
                disp_[vert][static_cast<std::size_t>(d)] *= heat_[vert];
            if(disp_norm_[vert] > 0.0){
                for(int d = 0; d < dim_; d++)
                    disp_[vert][static_cast<std::size_t>(d)] /= disp_norm_[vert];
            }
            if(misf_level == 0 &&
               final_stage_mode_ == FINAL_STAGE_FR_ND &&
               round_in_level > 1 &&
               final_move_scale_after_first_ < 1.0){
                for(int d = 0; d < dim_; d++)
                    disp_[vert][static_cast<std::size_t>(d)] *=
                        final_move_scale_after_first_;
            }
            if(record_steps){
                step_rows[static_cast<std::size_t>(i)] =
                    record_refinement_step_pre(vert,
                                               i + 1,
                                               active_count,
                                               level_index,
                                               misf_level,
                                               round_in_level,
                                               coord_before,
                                               pre_temp_disp,
                                               last_attraction_disp_,
                                               last_repulsion_disp_,
                                               disp_[vert],
                                               heat_before,
                                               heat_after,
                                               old_cos_before,
                                               old_cos_after,
                                               old_disp_norm_before,
                                               pre_temp_disp_norm,
                                               last_attraction_edges_,
                                               last_repulsion_neighbors_);
            }
        }

        for(int i = 0; i < active_count; i++){
            const vertex_t vert = misf.order[static_cast<std::size_t>(i)];
            for(int d = 0; d < dim_; d++){
                const std::size_t dd = static_cast<std::size_t>(d);
                coords_[vert][dd] += disp_[vert][dd];
                if(!std::isfinite(coords_[vert][dd]))
                    coords_[vert][dd] = 0.0;
            }
            const std::size_t row = step_rows[static_cast<std::size_t>(i)];
            if(record_steps && row != static_cast<std::size_t>(-1))
                record_refinement_step_after(row, vert);
        }

        if(trace != nullptr && ((round + 1) % trace_every) == 0)
            record_trace(trace, "round", round + 1, active_count, level_index, misf_level);
    }
}

std::vector<MetricNeighborND> DrawGraphND::metric_neighbors_for_level(
    vertex_t root,
    const WeightedMisfND &misf,
    int target_level) const
{
    if(target_level < 0 ||
       root >= misf.vertex_depth.size() ||
       misf.vertex_depth[root] < target_level)
        return std::vector<MetricNeighborND>();

    populate_metric_neighbors_for_root(root, misf);
    if(root >= metric_neighbors_cache_.size() ||
       target_level >= static_cast<int>(metric_neighbors_cache_[root].size()))
        return std::vector<MetricNeighborND>();
    return metric_neighbors_cache_[root][static_cast<std::size_t>(target_level)];
}

void DrawGraphND::populate_metric_neighbors_for_root(
    vertex_t root,
    const WeightedMisfND &misf) const
{
    if(root >= graph_.size() ||
       root >= misf.vertex_depth.size() ||
       root >= metric_neighbors_cache_.size())
        return;

    if(metric_neighbors_cached_[root])
        return;

    const int root_depth = misf.vertex_depth[root];
    metric_neighbors_cache_[root].assign(static_cast<std::size_t>(root_depth + 1),
                                         std::vector<MetricNeighborND>());
    metric_neighbors_cached_[root] = 1;
    if(root_depth < 0)
        return;

    const double inf = std::numeric_limits<double>::infinity();
    struct Node {
        double dist;
        vertex_t vert;
    };
    struct Greater {
        bool operator()(const Node &lhs, const Node &rhs) const
        {
            if(lhs.dist != rhs.dist)
                return lhs.dist > rhs.dist;
            return lhs.vert > rhs.vert;
        }
    };

    std::vector<double> dist(graph_.size(), inf);
    std::priority_queue<Node, std::vector<Node>, Greater> pq;
    dist[root] = 0.0;
    pq.push(Node{0.0, root});
    int bottom_level = 0;
    int settled = 0;
    const int max_settled = metric_neighbor_cap_for_root(root, misf);
    while(!pq.empty() && bottom_level <= root_depth){
        Node node = pq.top();
        pq.pop();
        if(node.dist > dist[node.vert] + kMetricTolND)
            continue;
        if(node.vert != root){
            settled++;
            if(node.vert < misf.vertex_depth.size() &&
               node.vert < graph_.size()){
                const int limit = std::min(misf.vertex_depth[node.vert], root_depth);
                for(int level = bottom_level; level <= limit; level++){
                    int target_count = num_nbrs_;
                    if(level < static_cast<int>(misf.level_size.size()))
                        target_count = std::min(
                            target_count,
                            std::max(0,
                                     misf.level_size[static_cast<std::size_t>(level)] - 1));
                    std::vector<MetricNeighborND> &layer =
                        metric_neighbors_cache_[root][static_cast<std::size_t>(level)];
                    if(static_cast<int>(layer.size()) < target_count){
                        layer.push_back(MetricNeighborND{node.vert, node.dist});
                    } else {
                        bottom_level = level + 1;
                    }
                }
            }
            if(max_settled > 0 && settled >= max_settled)
                break;
        }

        for(const NeighborND &nbr : graph_.neighbors(node.vert)){
            const double alt = node.dist + nbr.weight;
            const double best = dist[nbr.vertex];
            const double scale = std::max(1.0,
                                          std::max(std::fabs(alt),
                                                   std::isfinite(best) ? std::fabs(best) : 0.0));
            if(!std::isfinite(best) || alt + kMetricTolND * scale < best){
                dist[nbr.vertex] = alt;
                pq.push(Node{alt, nbr.vertex});
            }
        }
    }
}

void DrawGraphND::legacy_weighted_kk_displacement(
    vertex_t vert,
    const std::vector<MetricNeighborND> &neighbors,
    int active_count,
    int misf_level)
{
    disp_[vert].fill(0.0);
    std::fill(last_attraction_disp_.begin(), last_attraction_disp_.end(), 0.0);
    std::fill(last_repulsion_disp_.begin(), last_repulsion_disp_.end(), 0.0);
    last_attraction_edges_.clear();
    last_repulsion_neighbors_.clear();
    last_attraction_term_neighbors_.clear();
    last_attraction_term_weights_.clear();
    last_attraction_term_norm2_.clear();
    last_attraction_term_desired_.clear();
    last_attraction_term_desired2_.clear();
    last_attraction_term_scale_.clear();
    last_attraction_term_delta_.clear();
    last_attraction_term_step_.clear();
    last_attraction_term_cumulative_.clear();
    std::ostringstream attraction_edges;
    attraction_edges << std::setprecision(17);
    for(const MetricNeighborND &neighbor : neighbors){
        if(neighbor.vert >= coords_.size() || neighbor.vert == vert || neighbor.dist <= 0.0)
            continue;
        double norm2 = 0.0;
        PointND vect(static_cast<std::size_t>(dim_));
        std::vector<double> delta_flat(static_cast<std::size_t>(dim_));
        std::vector<double> step_flat(static_cast<std::size_t>(dim_));
        std::vector<double> cumulative_flat(static_cast<std::size_t>(dim_));
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] = coords_[neighbor.vert][dd] - coords_[vert][dd];
            norm2 += vect[dd] * vect[dd];
            delta_flat[dd] = vect[dd];
        }
        const double desired = neighbor.dist * kLegacyEdgeND;
        const double desired2 = desired * desired;
        const double scale = norm2 / desired2 - 1.0;
        if(attraction_edges.tellp() > 0)
            attraction_edges << ";";
        attraction_edges << static_cast<int>(neighbor.vert) + 1 << ":"
                         << neighbor.dist;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] *= scale;
            disp_[vert][dd] += vect[dd];
            last_attraction_disp_[dd] += vect[dd];
            step_flat[dd] = vect[dd];
            cumulative_flat[dd] = last_attraction_disp_[dd];
        }
        last_attraction_term_neighbors_.push_back(static_cast<int>(neighbor.vert) + 1);
        last_attraction_term_weights_.push_back(neighbor.dist);
        last_attraction_term_norm2_.push_back(norm2);
        last_attraction_term_desired_.push_back(desired);
        last_attraction_term_desired2_.push_back(desired2);
        last_attraction_term_scale_.push_back(scale);
        last_attraction_term_delta_.push_back(std::move(delta_flat));
        last_attraction_term_step_.push_back(std::move(step_flat));
        last_attraction_term_cumulative_.push_back(std::move(cumulative_flat));
    }
    last_attraction_edges_ = attraction_edges.str();

    if(misf_level > 0 && active_count > 1)
        add_legacy_active_repulsion(
            vert,
            active_count,
            coarse_repulsion_factor_ * 0.05 * kLegacyEdgeND * kLegacyEdgeND
        );

    scale_legacy_displacement(vert);
}

void DrawGraphND::legacy_weighted_kk_final_displacement(
    vertex_t vert,
    const std::vector<MetricNeighborND> &neighbors,
    int active_count)
{
    disp_[vert].fill(0.0);
    std::fill(last_attraction_disp_.begin(), last_attraction_disp_.end(), 0.0);
    std::fill(last_repulsion_disp_.begin(), last_repulsion_disp_.end(), 0.0);
    last_attraction_edges_.clear();
    last_repulsion_neighbors_.clear();
    last_attraction_term_neighbors_.clear();
    last_attraction_term_weights_.clear();
    last_attraction_term_norm2_.clear();
    last_attraction_term_desired_.clear();
    last_attraction_term_desired2_.clear();
    last_attraction_term_scale_.clear();
    last_attraction_term_delta_.clear();
    last_attraction_term_step_.clear();
    last_attraction_term_cumulative_.clear();
    std::ostringstream attraction_edges;
    attraction_edges << std::setprecision(17);

    for(const MetricNeighborND &neighbor : neighbors){
        if(neighbor.vert >= coords_.size() || neighbor.vert == vert || neighbor.dist <= 0.0)
            continue;
        double norm2 = 0.0;
        PointND vect(static_cast<std::size_t>(dim_));
        std::vector<double> delta_flat(static_cast<std::size_t>(dim_));
        std::vector<double> step_flat(static_cast<std::size_t>(dim_));
        std::vector<double> cumulative_flat(static_cast<std::size_t>(dim_));
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] = coords_[neighbor.vert][dd] - coords_[vert][dd];
            norm2 += vect[dd] * vect[dd];
            delta_flat[dd] = vect[dd];
        }
        const double desired = neighbor.dist * kLegacyEdgeND;
        const double desired2 = desired * desired;
        const double scale = norm2 / desired2 - 1.0;
        if(attraction_edges.tellp() > 0)
            attraction_edges << ";";
        attraction_edges << static_cast<int>(neighbor.vert) + 1 << ":"
                         << neighbor.dist;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] *= scale;
            disp_[vert][dd] += vect[dd];
            last_attraction_disp_[dd] += vect[dd];
            step_flat[dd] = vect[dd];
            cumulative_flat[dd] = last_attraction_disp_[dd];
        }
        last_attraction_term_neighbors_.push_back(static_cast<int>(neighbor.vert) + 1);
        last_attraction_term_weights_.push_back(neighbor.dist);
        last_attraction_term_norm2_.push_back(norm2);
        last_attraction_term_desired_.push_back(desired);
        last_attraction_term_desired2_.push_back(desired2);
        last_attraction_term_scale_.push_back(scale);
        last_attraction_term_delta_.push_back(std::move(delta_flat));
        last_attraction_term_step_.push_back(std::move(step_flat));
        last_attraction_term_cumulative_.push_back(std::move(cumulative_flat));
    }
    last_attraction_edges_ = attraction_edges.str();

    if(active_count > 1)
        add_legacy_active_repulsion(
            vert,
            active_count,
            repulsion_factor_ * 0.05 * kLegacyEdgeND * kLegacyEdgeND
        );

    scale_legacy_displacement(vert);
}

void DrawGraphND::legacy_weighted_fr_displacement(
    vertex_t vert,
    const std::vector<MetricNeighborND> &neighbors,
    int active_count)
{
    (void)active_count;
    disp_[vert].fill(0.0);
    std::fill(last_attraction_disp_.begin(), last_attraction_disp_.end(), 0.0);
    std::fill(last_repulsion_disp_.begin(), last_repulsion_disp_.end(), 0.0);
    last_attraction_term_neighbors_.clear();
    last_attraction_term_weights_.clear();
    last_attraction_term_norm2_.clear();
    last_attraction_term_desired_.clear();
    last_attraction_term_desired2_.clear();
    last_attraction_term_scale_.clear();
    last_attraction_term_delta_.clear();
    last_attraction_term_step_.clear();
    last_attraction_term_cumulative_.clear();
    std::ostringstream attraction_edges;
    std::ostringstream repulsion_neighbors;
    attraction_edges << std::setprecision(17);
    repulsion_neighbors << std::setprecision(17);

    for(const NeighborND &neighbor : graph_.neighbors(vert)){
        const vertex_t overt = neighbor.vertex;
        if(overt >= coords_.size())
            continue;
        double norm2 = 0.0;
        PointND vect(static_cast<std::size_t>(dim_));
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] = coords_[overt][dd] - coords_[vert][dd];
            norm2 += vect[dd] * vect[dd];
        }
        const double desired = kLegacyEdgeND * neighbor.weight;
        const double desired2 = desired * desired;
        if(desired2 <= 0.0)
            continue;
        const double scale = norm2 / desired2;
        std::vector<double> delta_flat(static_cast<std::size_t>(dim_));
        std::vector<double> step_flat(static_cast<std::size_t>(dim_));
        std::vector<double> cumulative_flat(static_cast<std::size_t>(dim_));
        if(attraction_edges.tellp() > 0)
            attraction_edges << ";";
        attraction_edges << static_cast<int>(overt) + 1 << ":" << neighbor.weight;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            delta_flat[dd] = vect[dd];
            vect[dd] *= norm2;
            vect[dd] /= desired2;
            disp_[vert][dd] += vect[dd];
            last_attraction_disp_[dd] += vect[dd];
            step_flat[dd] = vect[dd];
            cumulative_flat[dd] = last_attraction_disp_[dd];
        }
        last_attraction_term_neighbors_.push_back(static_cast<int>(overt) + 1);
        last_attraction_term_weights_.push_back(neighbor.weight);
        last_attraction_term_norm2_.push_back(norm2);
        last_attraction_term_desired_.push_back(desired);
        last_attraction_term_desired2_.push_back(desired2);
        last_attraction_term_scale_.push_back(scale);
        last_attraction_term_delta_.push_back(std::move(delta_flat));
        last_attraction_term_step_.push_back(std::move(step_flat));
        last_attraction_term_cumulative_.push_back(std::move(cumulative_flat));
    }

    const double fedge2 = repulsion_factor_ * 0.05 * kLegacyEdgeND * kLegacyEdgeND;
    for(const MetricNeighborND &neighbor : neighbors){
        if(neighbor.vert >= coords_.size() || neighbor.vert == vert || neighbor.dist <= 0.0)
            continue;
        double norm2 = 0.0;
        PointND vect(static_cast<std::size_t>(dim_));
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] = coords_[vert][dd] - coords_[neighbor.vert][dd];
            norm2 += vect[dd] * vect[dd];
        }
        if(norm2 <= 0.0)
            continue;
        const double scale = fedge2 / norm2;
        if(repulsion_neighbors.tellp() > 0)
            repulsion_neighbors << ";";
        repulsion_neighbors << static_cast<int>(neighbor.vert) + 1 << ":" << neighbor.dist;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            const double step = vect[dd] * scale;
            disp_[vert][dd] += step;
            last_repulsion_disp_[dd] += step;
        }
    }

    add_final_anchor_force(vert);
    scale_legacy_displacement(vert);
    last_attraction_edges_ = attraction_edges.str();
    last_repulsion_neighbors_ = repulsion_neighbors.str();
}

void DrawGraphND::add_legacy_active_repulsion(vertex_t vert,
                                              int active_count,
                                              double repulsion_scale)
{
    if(repulsion_scale <= 0.0 || active_count <= 1 ||
       coarse_repulsion_sample_ == 0)
        return;

    const int population = active_count - 1;
    if(active_count <= coarse_repulsion_exact_below_ ||
       coarse_repulsion_sample_ >= population){
        add_legacy_active_repulsion_exact(vert, active_count, repulsion_scale);
        return;
    }

    add_legacy_active_repulsion_sampled(
        vert,
        active_count,
        std::min(coarse_repulsion_sample_, population),
        repulsion_scale
    );
}

void DrawGraphND::add_legacy_active_repulsion_exact(vertex_t vert,
                                                    int active_count,
                                                    double repulsion_scale)
{
    std::ostringstream repulsion_neighbors;
    repulsion_neighbors << std::setprecision(17);
    if(!last_repulsion_neighbors_.empty())
        repulsion_neighbors << last_repulsion_neighbors_;
    for(int i = 0; i < active_count && i < static_cast<int>(trace_order_.size()); i++){
        const vertex_t overt = trace_order_[static_cast<std::size_t>(i)];
        if(overt == vert)
            continue;
        double norm2 = 0.0;
        PointND vect(static_cast<std::size_t>(dim_));
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] = coords_[vert][dd] - coords_[overt][dd];
            norm2 += vect[dd] * vect[dd];
        }
        if(norm2 <= 0.0)
            continue;
        const double scale = repulsion_scale / norm2;
        if(repulsion_neighbors.tellp() > 0)
            repulsion_neighbors << ";";
        repulsion_neighbors << static_cast<int>(overt) + 1;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] *= scale;
            disp_[vert][dd] += vect[dd];
            last_repulsion_disp_[dd] += vect[dd];
        }
    }
    last_repulsion_neighbors_ = repulsion_neighbors.str();
}

void DrawGraphND::add_legacy_active_repulsion_sampled(vertex_t vert,
                                                      int active_count,
                                                      int sample_count,
                                                      double repulsion_scale)
{
    if(sample_count <= 0)
        return;

    std::vector<vertex_t> sampled;
    sampled.reserve(static_cast<std::size_t>(sample_count));
    while(static_cast<int>(sampled.size()) < sample_count){
        const int idx = static_cast<int>(next_legacy_rand() %
                                         static_cast<unsigned long>(active_count));
        const vertex_t overt = trace_order_[static_cast<std::size_t>(idx)];
        if(overt == vert)
            continue;
        if(std::find(sampled.begin(), sampled.end(), overt) != sampled.end())
            continue;
        sampled.push_back(overt);
    }

    const double sample_scale =
        static_cast<double>(active_count - 1) / static_cast<double>(sample_count);
    std::ostringstream repulsion_neighbors;
    repulsion_neighbors << std::setprecision(17);
    if(!last_repulsion_neighbors_.empty())
        repulsion_neighbors << last_repulsion_neighbors_;
    for(vertex_t overt : sampled){
        double norm2 = 0.0;
        PointND vect(static_cast<std::size_t>(dim_));
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] = coords_[vert][dd] - coords_[overt][dd];
            norm2 += vect[dd] * vect[dd];
        }
        if(norm2 <= 0.0)
            continue;
        const double scale = (repulsion_scale * sample_scale) / norm2;
        if(repulsion_neighbors.tellp() > 0)
            repulsion_neighbors << ";";
        repulsion_neighbors << static_cast<int>(overt) + 1;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] *= scale;
            disp_[vert][dd] += vect[dd];
            last_repulsion_disp_[dd] += vect[dd];
        }
    }
    last_repulsion_neighbors_ = repulsion_neighbors.str();
}

void DrawGraphND::add_final_anchor_force(vertex_t vert)
{
    if(final_anchor_factor_ <= 0.0 || !final_anchor_ready_ ||
       vert >= final_anchor_pos_.size())
        return;

    for(int d = 0; d < dim_; d++){
        const std::size_t dd = static_cast<std::size_t>(d);
        disp_[vert][dd] +=
            (final_anchor_pos_[vert][dd] - coords_[vert][dd]) *
            (final_anchor_factor_ / kLegacyEdgeND);
    }
}

void DrawGraphND::prepare_final_anchors(int active_count)
{
    active_count = std::min(active_count, static_cast<int>(trace_order_.size()));
    for(int i = 0; i < active_count; i++){
        const vertex_t vert = trace_order_[static_cast<std::size_t>(i)];
        if(vert < final_anchor_pos_.size())
            final_anchor_pos_[vert] = coords_[vert];
    }
    final_anchor_ready_ = true;
}

void DrawGraphND::scale_legacy_displacement(vertex_t vert)
{
    const double norm = disp_[vert].norm();
    disp_norm_[vert] = round_legacy_norm_nd(norm);
    if(disp_norm_[vert] > 0.0 && norm > 0.0){
        const double scale = kLegacyEdgeND / norm;
        for(int d = 0; d < dim_; d++)
            disp_[vert][static_cast<std::size_t>(d)] *= scale;
        disp_norm_[vert] = disp_[vert].norm();
    }
}

void DrawGraphND::update_local_temperature(vertex_t vert)
{
    const double norm_old = old_disp_norm_[vert];
    const double norm_new = disp_norm_[vert];
    if(norm_old == 0.0 || norm_new == 0.0)
        return;

    const double cos = dot_point_nd(disp_[vert], old_disp_[vert]) / (norm_old * norm_new);
    if(old_cos_[vert] * cos > 0.0)
        heat_[vert] += heat_[vert] * s_ * cos * r_;
    else
        heat_[vert] += heat_[vert] * cos * r_;
    old_cos_[vert] = cos;
}

void DrawGraphND::initialize_multiscale_trace(const WeightedMisfND &misf)
{
    trace_order_ = misf.order;
    initialize_top_level(misf);
}

int DrawGraphND::lgkk_round_budget_for_layer(int misf_level) const
{
    int budget = 0;
    if(misf_level == 0)
        budget = lgkk_rounds_final_;
    else if(misf_level == 1)
        budget = lgkk_rounds_pre_final_;
    else
        budget = lgkk_rounds_coarse_;

    if(budget == 0)
        budget = lgkk_multiscale_rounds_;
    return budget;
}

bool DrawGraphND::should_run_multiscale_lgkk(int active_count,
                                             int misf_level) const
{
    if(lgkk_round_budget_for_layer(misf_level) == 0)
        return false;
    if(lgkk_local_nbrs_ == 0 && lgkk_landmark_count_ == 0)
        return false;
    if(active_count < 2 || active_count > lgkk_active_limit_)
        return false;
    if(lgkk_scope_ == LGKK_SCOPE_COARSE_ND && misf_level == 0)
        return false;
    return true;
}

void DrawGraphND::clear_lgkk_level_cache()
{
    lgkk_cache_active_count_ = 0;
    lgkk_cache_misf_level_ = -1;
    lgkk_cache_scale_l0_ = 1.0;
    lgkk_active_index_.clear();
    lgkk_distance_matrix_.clear();
    lgkk_pairs_.clear();
}

void DrawGraphND::compute_lgkk_active_shortest_paths(
    const WeightedMisfND &misf,
    int source_index,
    int active_count,
    std::vector<double> &dist,
    std::vector<int> *parent)
{
    const double inf = std::numeric_limits<double>::infinity();
    const double tol = 1e-10;
    dist.assign(static_cast<std::size_t>(active_count), inf);
    if(parent != nullptr)
        parent->assign(static_cast<std::size_t>(active_count), -1);
    if(source_index < 0 || source_index >= active_count)
        return;

    struct QueueNode {
        double dist;
        vertex_t vert;
        int index;
    };
    struct QueueNodeGreater {
        bool operator()(const QueueNode &lhs, const QueueNode &rhs) const
        {
            if(lhs.dist != rhs.dist)
                return lhs.dist > rhs.dist;
            return lhs.vert > rhs.vert;
        }
    };

    dist[static_cast<std::size_t>(source_index)] = 0.0;
    std::priority_queue<QueueNode,
                        std::vector<QueueNode>,
                        QueueNodeGreater> pq;
    pq.push(QueueNode{
        0.0,
        misf.order[static_cast<std::size_t>(source_index)],
        source_index
    });

    while(!pq.empty()){
        const double current_dist = pq.top().dist;
        const int current_index = pq.top().index;
        pq.pop();
        if(current_dist > dist[static_cast<std::size_t>(current_index)] + tol)
            continue;
        const vertex_t current_vert = misf.order[static_cast<std::size_t>(current_index)];
        for(const NeighborND &neighbor : graph_.neighbors(current_vert)){
            const vertex_t overt = neighbor.vertex;
            if(overt >= lgkk_active_index_.size())
                continue;
            const int overt_index = lgkk_active_index_[overt];
            if(overt_index < 0 || overt_index >= active_count)
                continue;
            const double alt = current_dist + neighbor.weight;
            const double best = dist[static_cast<std::size_t>(overt_index)];
            const double scale = std::max(1.0,
                                          std::max(std::fabs(alt),
                                                   std::isfinite(best) ? std::fabs(best) : 0.0));
            const bool improve =
                !std::isfinite(best) || alt + tol * scale < best;
            const bool equal =
                std::isfinite(best) && std::fabs(alt - best) <= tol * scale;
            if(improve){
                dist[static_cast<std::size_t>(overt_index)] = alt;
                if(parent != nullptr)
                    (*parent)[static_cast<std::size_t>(overt_index)] = current_index;
                pq.push(QueueNode{alt, overt, overt_index});
            } else if(equal && parent != nullptr){
                const int current_parent =
                    (*parent)[static_cast<std::size_t>(overt_index)];
                if(current_parent < 0 ||
                   misf.order[static_cast<std::size_t>(current_index)] <
                       misf.order[static_cast<std::size_t>(current_parent)]){
                    (*parent)[static_cast<std::size_t>(overt_index)] = current_index;
                }
            }
        }
    }
}

std::vector<int> DrawGraphND::lgkk_choose_local_neighbors(
    int source_index,
    int active_count) const
{
    struct Candidate {
        double dist;
        int index;
        vertex_t vert;
    };
    if(lgkk_local_nbrs_ == 0 || source_index < 0 || source_index >= active_count)
        return std::vector<int>();

    const double *row =
        lgkk_distance_matrix_.data() +
        static_cast<std::size_t>(source_index * active_count);
    std::vector<Candidate> candidates;
    candidates.reserve(static_cast<std::size_t>(std::max(0, active_count - 1)));
    for(int idx = 0; idx < active_count; idx++){
        if(idx == source_index || !std::isfinite(row[idx]))
            continue;
        candidates.push_back(Candidate{
            row[idx],
            idx,
            trace_order_[static_cast<std::size_t>(idx)]
        });
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate &lhs, const Candidate &rhs){
                  if(lhs.dist != rhs.dist)
                      return lhs.dist < rhs.dist;
                  return lhs.vert < rhs.vert;
              });

    std::vector<int> out;
    out.reserve(static_cast<std::size_t>(
        std::min<int>(lgkk_local_nbrs_, static_cast<int>(candidates.size()))));
    for(std::size_t i = 0;
        i < candidates.size() && static_cast<int>(out.size()) < lgkk_local_nbrs_;
        i++){
        out.push_back(candidates[i].index);
    }
    return out;
}

std::vector<int> DrawGraphND::lgkk_choose_landmarks(
    const WeightedMisfND &misf,
    int source_index,
    int active_count) const
{
    if(lgkk_landmark_count_ == 0 || source_index < 0 || source_index >= active_count)
        return std::vector<int>();

    const double *source_row =
        lgkk_distance_matrix_.data() +
        static_cast<std::size_t>(source_index * active_count);
    std::vector<int> candidates;
    candidates.reserve(static_cast<std::size_t>(std::max(0, active_count - 1)));
    for(int idx = 0; idx < active_count; idx++){
        if(idx == source_index || !std::isfinite(source_row[idx]))
            continue;
        candidates.push_back(idx);
    }
    if(candidates.empty())
        return std::vector<int>();

    std::vector<int> selected;
    std::vector<double> coverage;
    coverage.reserve(candidates.size());
    for(std::size_t i = 0; i < candidates.size(); i++)
        coverage.push_back(source_row[candidates[i]]);

    const std::size_t max_count =
        std::min<std::size_t>(static_cast<std::size_t>(lgkk_landmark_count_),
                              candidates.size());
    for(std::size_t step = 0; step < max_count; step++){
        std::size_t choice_pos = candidates.size();
        double choice_score = -1.0;
        for(std::size_t i = 0; i < candidates.size(); i++){
            const double score = (step == 0) ? source_row[candidates[i]] : coverage[i];
            if(score > choice_score ||
               (std::fabs(score - choice_score) <= 1e-12 &&
                (choice_pos >= candidates.size() ||
                 misf.order[static_cast<std::size_t>(candidates[i])] <
                     misf.order[static_cast<std::size_t>(candidates[choice_pos])]))){
                choice_score = score;
                choice_pos = i;
            }
        }
        if(choice_pos >= candidates.size())
            break;
        const int choice = candidates[choice_pos];
        selected.push_back(choice);
        candidates.erase(candidates.begin() + static_cast<long>(choice_pos));
        coverage.erase(coverage.begin() + static_cast<long>(choice_pos));
        if(candidates.empty())
            break;

        const double *choice_row =
            lgkk_distance_matrix_.data() +
            static_cast<std::size_t>(choice * active_count);
        for(std::size_t i = 0; i < candidates.size(); i++)
            coverage[i] = std::min(coverage[i], choice_row[candidates[i]]);
    }
    return selected;
}

void DrawGraphND::build_lgkk_level_cache(const WeightedMisfND &misf,
                                         int active_count,
                                         int misf_level)
{
    clear_lgkk_level_cache();
    if(!should_run_multiscale_lgkk(active_count, misf_level))
        return;

    active_count = std::min(active_count, static_cast<int>(misf.order.size()));
    lgkk_active_index_.assign(graph_.size(), -1);
    for(int i = 0; i < active_count; i++)
        lgkk_active_index_[misf.order[static_cast<std::size_t>(i)]] = i;

    lgkk_distance_matrix_.assign(
        static_cast<std::size_t>(active_count * active_count),
        std::numeric_limits<double>::infinity());
    std::vector<double> dist;
    for(int source_index = 0; source_index < active_count; source_index++){
        compute_lgkk_active_shortest_paths(misf,
                                           source_index,
                                           active_count,
                                           dist,
                                           nullptr);
        std::copy(dist.begin(),
                  dist.end(),
                  lgkk_distance_matrix_.begin() +
                      static_cast<std::ptrdiff_t>(source_index * active_count));
    }

    std::vector<std::vector<int>> selected_targets(static_cast<std::size_t>(active_count));
    for(int source_index = 0; source_index < active_count; source_index++){
        std::vector<int> local =
            lgkk_choose_local_neighbors(source_index, active_count);
        std::vector<int> landmarks =
            lgkk_choose_landmarks(misf, source_index, active_count);
        local.insert(local.end(), landmarks.begin(), landmarks.end());
        std::sort(local.begin(), local.end());
        local.erase(std::unique(local.begin(), local.end()), local.end());
        selected_targets[static_cast<std::size_t>(source_index)] = std::move(local);
    }

    std::unordered_set<std::uint64_t> seen_pairs;
    std::vector<int> parent;
    lgkk_pairs_.reserve(
        static_cast<std::size_t>(active_count) *
        static_cast<std::size_t>(std::max(1, lgkk_local_nbrs_ + lgkk_landmark_count_)));
    for(int source_index = 0; source_index < active_count; source_index++){
        if(selected_targets[static_cast<std::size_t>(source_index)].empty())
            continue;
        compute_lgkk_active_shortest_paths(misf,
                                           source_index,
                                           active_count,
                                           dist,
                                           &parent);
        for(int target_index : selected_targets[static_cast<std::size_t>(source_index)]){
            const vertex_t source_vert =
                misf.order[static_cast<std::size_t>(source_index)];
            const vertex_t target_vert =
                misf.order[static_cast<std::size_t>(target_index)];
            const vertex_t min_vert = std::min(source_vert, target_vert);
            const vertex_t max_vert = std::max(source_vert, target_vert);
            const std::uint64_t key =
                (static_cast<std::uint64_t>(min_vert) << 32) |
                static_cast<std::uint64_t>(max_vert);
            if(seen_pairs.find(key) != seen_pairs.end())
                continue;
            if(target_index < 0 ||
               target_index >= static_cast<int>(dist.size()) ||
               !std::isfinite(dist[static_cast<std::size_t>(target_index)]))
                continue;

            std::vector<vertex_t> path_vertices;
            int current_index = target_index;
            while(current_index >= 0 && current_index != source_index){
                path_vertices.push_back(
                    misf.order[static_cast<std::size_t>(current_index)]);
                current_index = parent[static_cast<std::size_t>(current_index)];
            }
            if(current_index < 0)
                continue;
            path_vertices.push_back(source_vert);
            std::reverse(path_vertices.begin(), path_vertices.end());
            if(path_vertices.size() < 2)
                continue;

            LgkkPairCacheND pair;
            pair.source = source_vert;
            pair.target = target_vert;
            pair.graph_distance = dist[static_cast<std::size_t>(target_index)];
            pair.path_edges.reserve(path_vertices.size() - 1);
            for(std::size_t edge_index = 1;
                edge_index < path_vertices.size();
                edge_index++){
                pair.path_edges.push_back(LgkkPathEdgeND{
                    path_vertices[edge_index - 1],
                    path_vertices[edge_index]
                });
            }
            lgkk_pairs_.push_back(std::move(pair));
            seen_pairs.insert(key);
        }
    }

    double numerator = 0.0;
    double denominator = 0.0;
    const double eps2 = 1e-16;
    for(const LgkkPairCacheND &pair : lgkk_pairs_){
        double h = 0.0;
        for(const LgkkPathEdgeND &edge_ref : pair.path_edges){
            double norm2 = 0.0;
            for(int d = 0; d < dim_; d++){
                const std::size_t dd = static_cast<std::size_t>(d);
                const double delta = coords_[edge_ref.u][dd] - coords_[edge_ref.v][dd];
                norm2 += delta * delta;
            }
            h += std::sqrt(norm2 + eps2);
        }
        const double g = std::max(pair.graph_distance, 1e-8);
        const double kk = 1.0 / (g * g);
        numerator += kk * g * h;
        denominator += kk * g * g;
    }

    lgkk_cache_scale_l0_ = denominator > 0.0 ? numerator / denominator : 1.0;
    lgkk_cache_active_count_ = active_count;
    lgkk_cache_misf_level_ = misf_level;
}

int DrawGraphND::lgkk_refine_level(const WeightedMisfND &misf,
                                   int active_count,
                                   int misf_level,
                                   int base_rounds,
                                   LayoutTraceND *trace,
                                   int trace_every,
                                   int level_index)
{
    if(!should_run_multiscale_lgkk(active_count, misf_level))
        return base_rounds;
    const int round_budget = lgkk_round_budget_for_layer(misf_level);
    if(round_budget == 0)
        return base_rounds;
    build_lgkk_level_cache(misf, active_count, misf_level);
    if(lgkk_pairs_.empty())
        return base_rounds;

    const double eps2 = 1e-16;
    const double initial_step = 1.0;
    const double step_shrink = 0.5;
    const double armijo = 1e-4;
    const double grad_tol2 = 1e-16;
    const double min_step = 1e-8;
    const double distance_floor = 1e-8;

    struct LgkkStateND {
        double energy;
        double grad_norm2;
        std::vector<PointND> gradient;
    };

    const std::size_t dim = static_cast<std::size_t>(dim_);
    active_count = std::min(active_count, static_cast<int>(misf.order.size()));
    std::vector<PointND> active_pos(static_cast<std::size_t>(active_count),
                                    PointND(dim));
    for(int i = 0; i < active_count; i++)
        active_pos[static_cast<std::size_t>(i)] =
            coords_[misf.order[static_cast<std::size_t>(i)]];

    auto evaluate_state = [&](const std::vector<PointND> &coords){
        LgkkStateND state;
        state.energy = 0.0;
        state.grad_norm2 = 0.0;
        state.gradient.assign(static_cast<std::size_t>(active_count), PointND(dim));
        for(PointND &grad : state.gradient)
            grad.fill(0.0);

        for(const LgkkPairCacheND &pair : lgkk_pairs_){
            const double g = std::max(pair.graph_distance, distance_floor);
            const double kk = 1.0 / (g * g);
            const double target = lgkk_cache_scale_l0_ * pair.graph_distance;

            std::vector<PointND> edge_diffs;
            std::vector<double> edge_lens;
            edge_diffs.reserve(pair.path_edges.size());
            edge_lens.reserve(pair.path_edges.size());

            double h = 0.0;
            for(const LgkkPathEdgeND &edge_ref : pair.path_edges){
                const int u_index = lgkk_active_index_[edge_ref.u];
                const int v_index = lgkk_active_index_[edge_ref.v];
                if(u_index < 0 || v_index < 0)
                    continue;
                PointND diff(dim);
                double norm2 = 0.0;
                for(int d = 0; d < dim_; d++){
                    const std::size_t dd = static_cast<std::size_t>(d);
                    diff[dd] =
                        coords[static_cast<std::size_t>(u_index)][dd] -
                        coords[static_cast<std::size_t>(v_index)][dd];
                    norm2 += diff[dd] * diff[dd];
                }
                const double len = std::sqrt(norm2 + eps2);
                edge_diffs.push_back(diff);
                edge_lens.push_back(len);
                h += len;
            }
            if(edge_diffs.empty())
                continue;

            const double resid = h - target;
            const double coeff = kk * resid;
            state.energy += 0.5 * kk * resid * resid;
            for(std::size_t edge_index = 0;
                edge_index < pair.path_edges.size();
                edge_index++){
                const LgkkPathEdgeND &edge_ref = pair.path_edges[edge_index];
                const int u_index = lgkk_active_index_[edge_ref.u];
                const int v_index = lgkk_active_index_[edge_ref.v];
                if(u_index < 0 || v_index < 0)
                    continue;
                if(edge_lens[edge_index] <= 0.0)
                    continue;
                const double scale = coeff / edge_lens[edge_index];
                for(int d = 0; d < dim_; d++){
                    const std::size_t dd = static_cast<std::size_t>(d);
                    const double step_vec = edge_diffs[edge_index][dd] * scale;
                    state.gradient[static_cast<std::size_t>(u_index)][dd] += step_vec;
                    state.gradient[static_cast<std::size_t>(v_index)][dd] -= step_vec;
                }
            }
        }

        for(const PointND &grad : state.gradient)
            state.grad_norm2 += grad.norm2();
        return state;
    };

    std::vector<PointND> accepted_move(static_cast<std::size_t>(active_count),
                                       PointND(dim));
    for(PointND &move : accepted_move)
        move.fill(0.0);

    LgkkStateND state = evaluate_state(active_pos);
    int trace_round_in_level = base_rounds;
    for(int round_index = 1; round_index <= round_budget; round_index++){
        if(!std::isfinite(state.energy) || state.grad_norm2 <= grad_tol2)
            break;

        double step = initial_step;
        bool accepted = false;
        std::vector<PointND> proposal(static_cast<std::size_t>(active_count),
                                      PointND(dim));
        LgkkStateND candidate = state;

        while(std::isfinite(step) && step >= min_step){
            for(int i = 0; i < active_count; i++){
                for(int d = 0; d < dim_; d++){
                    const std::size_t dd = static_cast<std::size_t>(d);
                    proposal[static_cast<std::size_t>(i)][dd] =
                        active_pos[static_cast<std::size_t>(i)][dd] -
                        state.gradient[static_cast<std::size_t>(i)][dd] * step;
                }
            }

            candidate = evaluate_state(proposal);
            const double target_energy =
                state.energy - armijo * step * state.grad_norm2;
            if(std::isfinite(candidate.energy) &&
               candidate.energy <= target_energy){
                accepted = true;
                break;
            }
            step *= step_shrink;
        }

        if(!accepted)
            break;

        for(int i = 0; i < active_count; i++){
            const std::size_t ii = static_cast<std::size_t>(i);
            const vertex_t vert = misf.order[ii];
            for(int d = 0; d < dim_; d++){
                const std::size_t dd = static_cast<std::size_t>(d);
                accepted_move[ii][dd] = proposal[ii][dd] - active_pos[ii][dd];
                active_pos[ii][dd] = proposal[ii][dd];
                coords_[vert][dd] = active_pos[ii][dd];
                disp_[vert][dd] = accepted_move[ii][dd];
                old_disp_[vert][dd] = accepted_move[ii][dd];
            }
            disp_norm_[vert] = round_legacy_norm_nd(accepted_move[ii].norm());
            old_disp_norm_[vert] = disp_norm_[vert];
        }

        state = candidate;
        trace_round_in_level = base_rounds + round_index;
        if(trace != nullptr && ((round_index % std::max(1, trace_every)) == 0))
            record_trace(trace,
                         "lgkk",
                         trace_round_in_level,
                         active_count,
                         level_index,
                         misf_level);
    }
    return trace_round_in_level;
}

void DrawGraphND::initialize_top_level(const WeightedMisfND &misf)
{
    for(PointND &point : coords_)
        point.fill(0.0);

    double diam = 0.0;
    for(int i = 0;
        i < misf.num_init && i < static_cast<int>(misf.order.size());
        i++){
        const vertex_t vert = misf.order[static_cast<std::size_t>(i)];
        diam = std::max(diam, std::ceil(weighted_eccentricity(vert)));
    }
    if(diam <= 0.0)
        diam = 1.0;

    box_size_ = kLegacyEdgeND * 0.7 * diam;
    box2_size_ = std::max(1, static_cast<int>(2.0 * box_size_ + 1.0));
    legacy_rng_state_ = misf.rng_state;
    PointND center(static_cast<std::size_t>(dim_));

    for(int i = 0; i < misf.num_init && i < static_cast<int>(misf.order.size()); i++){
        const vertex_t vert = misf.order[static_cast<std::size_t>(i)];
        for(int d = 0; d < dim_; d++){
            coords_[vert][static_cast<std::size_t>(d)] =
                static_cast<double>(next_legacy_rand() %
                                    static_cast<unsigned long>(box2_size_)) - box_size_;
            center[static_cast<std::size_t>(d)] += coords_[vert][static_cast<std::size_t>(d)];
        }
    }

    for(int d = 0; d < dim_; d++)
        center[static_cast<std::size_t>(d)] /= static_cast<double>(std::max(1, misf.num_init));
    for(int i = 0; i < misf.num_init && i < static_cast<int>(misf.order.size()); i++){
        const vertex_t vert = misf.order[static_cast<std::size_t>(i)];
        for(int d = 0; d < dim_; d++)
            coords_[vert][static_cast<std::size_t>(d)] -= center[static_cast<std::size_t>(d)];
    }
}

void DrawGraphND::insert_level_vertices(const WeightedMisfND &misf,
                                        int previous_active_count,
                                        int active_count,
                                        int level_index,
                                        int misf_level)
{
    previous_active_count = std::max(0, previous_active_count);
    active_count = std::min(active_count, static_cast<int>(misf.order.size()));
    for(int i = previous_active_count; i < active_count; i++){
        const vertex_t root = misf.order[static_cast<std::size_t>(i)];
        const int root_depth = misf.vertex_depth[root];
        const bool level0_insertion = root_depth == 0;
        const int anchor_count = level0_insertion ? level0_anchor_count_ : insertion_anchor_count_;
        const int local_kk_steps = level0_insertion ? level0_local_kk_steps_ : 3;
        int insertion_mode = level0_insertion ? level0_insertion_mode_ : LEVEL0_INSERT_INHERIT_ND;
        if(insertion_mode == LEVEL0_INSERT_INHERIT_ND)
            insertion_mode = placement_mode_;

        const std::vector<std::pair<vertex_t, double>> anchors =
            weighted_anchor_distances(root, misf, previous_active_count,
                                      root_depth, anchor_count,
                                      insertion_anchor_strategy_);
        std::vector<vertex_t> anchor_vertices;
        anchor_vertices.reserve(anchors.size());
        for(const auto &anchor : anchors)
            anchor_vertices.push_back(anchor.first);
        PointND coord_initial(static_cast<std::size_t>(dim_));
        PointND old_disp_initial(static_cast<std::size_t>(dim_));
        double old_disp_norm_initial = 0.0;
        if(anchor_vertices.empty()){
            coords_[root] = random_legacy_point();
            old_disp_[root].fill(0.0);
            old_disp_norm_[root] = 0.0;
            coord_initial = coords_[root];
            old_disp_initial = old_disp_[root];
            old_disp_norm_initial = old_disp_norm_[root];
        } else {
            if(insertion_mode == LEVEL0_INSERT_LEAST_SQUARES_ND)
                coords_[root] = weighted_least_squares_position(anchors);
            else if(insertion_mode == PLACEMENT_CIRCLE_ND && dim_ == 2)
                coords_[root] = weighted_circle_position(anchors);
            else
                coords_[root] = weighted_barycenter(anchors);
            coord_initial = coords_[root];

            old_disp_[root].fill(0.0);
            for(vertex_t anchor : anchor_vertices){
                for(int d = 0; d < dim_; d++)
                    old_disp_[root][static_cast<std::size_t>(d)] +=
                        old_disp_[anchor][static_cast<std::size_t>(d)];
            }
            for(int d = 0; d < dim_; d++)
                old_disp_[root][static_cast<std::size_t>(d)] /=
                    static_cast<double>(anchor_vertices.size());
            old_disp_norm_[root] = old_disp_[root].norm();
            old_disp_initial = old_disp_[root];
            old_disp_norm_initial = old_disp_norm_[root];

            polish_inserted_vertex(root, anchors, local_kk_steps);
        }
        record_insertion_trace(root,
                               level_index,
                               misf_level,
                               previous_active_count,
                               active_count,
                               i + 1,
                               root_depth,
                               anchor_count,
                               static_cast<int>(anchor_vertices.size()),
                               insertion_mode,
                               local_kk_steps,
                               anchors,
                               coord_initial,
                               coords_[root],
                               old_disp_initial,
                               old_disp_[root],
                               old_disp_norm_initial,
                               old_disp_norm_[root]);
    }
}

std::vector<std::pair<vertex_t, double>> DrawGraphND::select_insertion_anchor_subset(
    const std::vector<std::pair<vertex_t, double>> &anchors,
    int target_count) const
{
    if(target_count <= 0 || static_cast<int>(anchors.size()) <= target_count)
        return anchors;

    PointND candidate_centroid(static_cast<std::size_t>(dim_));
    for(const auto &anchor : anchors){
        for(int d = 0; d < dim_; d++)
            candidate_centroid[static_cast<std::size_t>(d)] +=
                coords_[anchor.first][static_cast<std::size_t>(d)];
    }
    for(int d = 0; d < dim_; d++)
        candidate_centroid[static_cast<std::size_t>(d)] /= static_cast<double>(anchors.size());

    auto point_dist2 = [&](vertex_t lhs, vertex_t rhs) {
        return squared_distance(coords_[lhs], coords_[rhs]);
    };
    auto centroid_dist2 = [&](vertex_t vert) {
        double out = 0.0;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            const double delta = coords_[vert][dd] - candidate_centroid[dd];
            out += delta * delta;
        }
        return out;
    };
    auto build_subset = [&](const std::vector<std::size_t> &indices) {
        std::vector<std::pair<vertex_t, double>> selected;
        selected.reserve(indices.size());
        for(std::size_t idx : indices)
            selected.push_back(anchors[idx]);
        return selected;
    };

    if(insertion_anchor_strategy_ == INSERT_ANCHOR_STRATEGY_SPREAD_PREV_ND &&
       dim_ == 2){
        const double pi = 3.14159265358979323846;
        struct PolarCandidate {
            std::size_t idx;
            double angle;
            double radius2;
            vertex_t vert;
        };
        auto angular_gap = [pi](double a, double b) {
            double diff = std::fabs(a - b);
            const double two_pi = 2.0 * pi;
            while(diff > two_pi)
                diff -= two_pi;
            if(diff > pi)
                diff = two_pi - diff;
            return diff;
        };

        std::vector<PolarCandidate> polar;
        polar.reserve(anchors.size());
        for(std::size_t i = 0; i < anchors.size(); i++){
            const PointND &p = coords_[anchors[i].first];
            const double dx = p[0] - candidate_centroid[0];
            const double dy = p[1] - candidate_centroid[1];
            polar.push_back(PolarCandidate{
                i,
                std::atan2(dy, dx),
                dx * dx + dy * dy,
                anchors[i].first
            });
        }

        std::vector<std::size_t> selected;
        std::vector<char> used(anchors.size(), 0);
        std::size_t first_choice = 0;
        for(std::size_t i = 1; i < polar.size(); i++){
            if(polar[i].radius2 > polar[first_choice].radius2 + 1e-12 ||
               (std::fabs(polar[i].radius2 - polar[first_choice].radius2) <= 1e-12 &&
                polar[i].vert < polar[first_choice].vert)){
                first_choice = i;
            }
        }
        selected.push_back(polar[first_choice].idx);
        used[polar[first_choice].idx] = 1;

        while(static_cast<int>(selected.size()) < target_count){
            std::size_t best_polar_index = polar.size();
            double best_angle_score = -1.0;
            double best_sep_score = -1.0;
            for(std::size_t i = 0; i < polar.size(); i++){
                if(used[polar[i].idx])
                    continue;
                double min_angle_gap = std::numeric_limits<double>::infinity();
                double min_sep = std::numeric_limits<double>::infinity();
                for(std::size_t chosen_idx : selected){
                    const PointND &chosen = coords_[anchors[chosen_idx].first];
                    const double cdx = chosen[0] - candidate_centroid[0];
                    const double cdy = chosen[1] - candidate_centroid[1];
                    const double chosen_angle = std::atan2(cdy, cdx);
                    min_angle_gap = std::min(min_angle_gap,
                                             angular_gap(polar[i].angle, chosen_angle));
                    min_sep = std::min(min_sep,
                                       point_dist2(anchors[polar[i].idx].first,
                                                   anchors[chosen_idx].first));
                }
                if(min_angle_gap > best_angle_score + 1e-12 ||
                   (std::fabs(min_angle_gap - best_angle_score) <= 1e-12 &&
                    (min_sep > best_sep_score + 1e-12 ||
                     (std::fabs(min_sep - best_sep_score) <= 1e-12 &&
                      (best_polar_index >= polar.size() ||
                       polar[i].vert < polar[best_polar_index].vert))))){
                    best_angle_score = min_angle_gap;
                    best_sep_score = min_sep;
                    best_polar_index = i;
                }
            }
            if(best_polar_index >= polar.size())
                break;
            used[polar[best_polar_index].idx] = 1;
            selected.push_back(polar[best_polar_index].idx);
        }

        if(static_cast<int>(selected.size()) == target_count)
            return build_subset(selected);
    }

    if(insertion_anchor_strategy_ == INSERT_ANCHOR_STRATEGY_BALANCED_BAND_ND){
        auto subset_objective = [&](const std::vector<std::size_t> &subset_indices) {
            PointND subset_centroid(static_cast<std::size_t>(dim_));
            for(std::size_t idx : subset_indices){
                for(int d = 0; d < dim_; d++)
                    subset_centroid[static_cast<std::size_t>(d)] +=
                        coords_[anchors[idx].first][static_cast<std::size_t>(d)];
            }
            for(int d = 0; d < dim_; d++)
                subset_centroid[static_cast<std::size_t>(d)] /=
                    static_cast<double>(subset_indices.size());

            double centroid_error = 0.0;
            for(int d = 0; d < dim_; d++){
                const std::size_t dd = static_cast<std::size_t>(d);
                const double delta = subset_centroid[dd] - candidate_centroid[dd];
                centroid_error += delta * delta;
            }

            double min_sep = 0.0;
            if(subset_indices.size() >= 2){
                min_sep = std::numeric_limits<double>::infinity();
                for(std::size_t i = 0; i < subset_indices.size(); i++){
                    for(std::size_t j = i + 1; j < subset_indices.size(); j++){
                        min_sep = std::min(min_sep,
                                           point_dist2(anchors[subset_indices[i]].first,
                                                       anchors[subset_indices[j]].first));
                    }
                }
                if(!std::isfinite(min_sep))
                    min_sep = 0.0;
            }

            double scale = 0.0;
            for(const auto &anchor : anchors)
                scale = std::max(scale, centroid_dist2(anchor.first));
            if(scale <= 0.0)
                scale = 1.0;
            return centroid_error / scale - 0.25 * (min_sep / scale);
        };

        const std::size_t n = anchors.size();
        const std::size_t k = static_cast<std::size_t>(std::min<int>(target_count, n));
        double combination_count = 1.0;
        for(std::size_t i = 1; i <= k; i++){
            combination_count *= static_cast<double>(n - k + i);
            combination_count /= static_cast<double>(i);
        }

        if(combination_count <= 50000.0){
            std::vector<std::size_t> current;
            std::vector<std::size_t> best;
            double best_objective = std::numeric_limits<double>::infinity();
            std::function<void(std::size_t, std::size_t)> dfs =
                [&](std::size_t start, std::size_t remaining) {
                    if(remaining == 0){
                        const double objective = subset_objective(current);
                        if(objective < best_objective - 1e-12){
                            best_objective = objective;
                            best = current;
                        } else if(std::fabs(objective - best_objective) <= 1e-12){
                            if(best.empty() || anchors[current[0]].first < anchors[best[0]].first)
                                best = current;
                        }
                        return;
                    }
                    const std::size_t limit = n - remaining;
                    for(std::size_t idx = start; idx <= limit; idx++){
                        current.push_back(idx);
                        dfs(idx + 1, remaining - 1);
                        current.pop_back();
                    }
                };
            dfs(0, k);
            if(!best.empty())
                return build_subset(best);
        }

        std::vector<std::size_t> greedy;
        std::vector<char> used(anchors.size(), 0);
        while(greedy.size() < k){
            std::size_t best_index = anchors.size();
            double best_objective = std::numeric_limits<double>::infinity();
            for(std::size_t i = 0; i < anchors.size(); i++){
                if(used[i])
                    continue;
                greedy.push_back(i);
                const double objective = subset_objective(greedy);
                greedy.pop_back();
                if(objective < best_objective - 1e-12 ||
                   (std::fabs(objective - best_objective) <= 1e-12 &&
                    (best_index >= anchors.size() ||
                     anchors[i].first < anchors[best_index].first))){
                    best_objective = objective;
                    best_index = i;
                }
            }
            if(best_index >= anchors.size())
                break;
            used[best_index] = 1;
            greedy.push_back(best_index);
        }
        if(!greedy.empty())
            return build_subset(greedy);
    }

    std::vector<std::size_t> selected;
    std::vector<char> used(anchors.size(), 0);
    std::size_t first_choice = 0;
    double first_score = -1.0;
    for(std::size_t i = 0; i < anchors.size(); i++){
        const double score = centroid_dist2(anchors[i].first);
        if(score > first_score ||
           (std::fabs(score - first_score) <= 1e-12 &&
            anchors[i].first < anchors[first_choice].first)){
            first_score = score;
            first_choice = i;
        }
    }
    selected.push_back(first_choice);
    used[first_choice] = 1;

    while(static_cast<int>(selected.size()) < target_count){
        std::size_t best_index = anchors.size();
        double best_score = -1.0;
        for(std::size_t i = 0; i < anchors.size(); i++){
            if(used[i])
                continue;
            double min_sep = std::numeric_limits<double>::infinity();
            for(std::size_t chosen_idx : selected)
                min_sep = std::min(
                    min_sep,
                    point_dist2(anchors[i].first, anchors[chosen_idx].first)
                );
            if(min_sep > best_score ||
               (std::fabs(min_sep - best_score) <= 1e-12 &&
                (best_index >= anchors.size() ||
                 anchors[i].first < anchors[best_index].first))){
                best_score = min_sep;
                best_index = i;
            }
        }
        if(best_index >= anchors.size())
            break;
        used[best_index] = 1;
        selected.push_back(best_index);
    }
    return build_subset(selected);
}

PointND DrawGraphND::weighted_barycenter(const std::vector<vertex_t> &anchors) const
{
    PointND out(static_cast<std::size_t>(dim_));
    if(anchors.empty())
        return out;
    for(int d = 0; d < dim_; d++)
        out[static_cast<std::size_t>(d)] =
            coords_[anchors[0]][static_cast<std::size_t>(d)];
    for(std::size_t i = 1; i < anchors.size(); i++){
        const vertex_t anchor = anchors[i];
        for(int d = 0; d < dim_; d++)
            out[static_cast<std::size_t>(d)] += coords_[anchor][static_cast<std::size_t>(d)];
    }
    for(int d = 0; d < dim_; d++)
        out[static_cast<std::size_t>(d)] /= static_cast<double>(anchors.size());
    return out;
}

PointND DrawGraphND::weighted_barycenter(
    const std::vector<std::pair<vertex_t, double>> &anchors) const
{
    std::vector<vertex_t> vertices;
    vertices.reserve(anchors.size());
    for(const auto &anchor : anchors)
        vertices.push_back(anchor.first);
    return weighted_barycenter(vertices);
}

PointND DrawGraphND::weighted_least_squares_position(
    const std::vector<std::pair<vertex_t, double>> &anchors) const
{
    if(static_cast<int>(anchors.size()) <= dim_)
        return weighted_barycenter(anchors);

    const PointND &base = coords_[anchors[0].first];
    const double d0 = anchors[0].second * kLegacyEdgeND;
    double base_norm2 = 0.0;
    for(int k = 0; k < dim_; k++)
        base_norm2 += base[static_cast<std::size_t>(k)] * base[static_cast<std::size_t>(k)];

    std::vector<std::vector<double>> ata(static_cast<std::size_t>(dim_),
                                         std::vector<double>(static_cast<std::size_t>(dim_), 0.0));
    std::vector<double> atb(static_cast<std::size_t>(dim_), 0.0);

    for(std::size_t i = 1; i < anchors.size(); i++){
        const PointND &anchor = coords_[anchors[i].first];
        const double di = anchors[i].second * kLegacyEdgeND;
        double anchor_norm2 = 0.0;
        std::vector<double> row(static_cast<std::size_t>(dim_), 0.0);
        for(int k = 0; k < dim_; k++){
            const std::size_t kk = static_cast<std::size_t>(k);
            row[kk] = 2.0 * (anchor[kk] - base[kk]);
            anchor_norm2 += anchor[kk] * anchor[kk];
        }
        const double rhs = anchor_norm2 - base_norm2 - di * di + d0 * d0;
        for(int r = 0; r < dim_; r++){
            const std::size_t rr = static_cast<std::size_t>(r);
            atb[rr] += row[rr] * rhs;
            for(int c = 0; c < dim_; c++){
                const std::size_t cc = static_cast<std::size_t>(c);
                ata[rr][cc] += row[rr] * row[cc];
            }
        }
    }

    std::vector<std::vector<double>> aug(
        static_cast<std::size_t>(dim_),
        std::vector<double>(static_cast<std::size_t>(dim_ + 1), 0.0)
    );
    for(int r = 0; r < dim_; r++){
        const std::size_t rr = static_cast<std::size_t>(r);
        for(int c = 0; c < dim_; c++)
            aug[rr][static_cast<std::size_t>(c)] = ata[rr][static_cast<std::size_t>(c)];
        aug[rr][static_cast<std::size_t>(dim_)] = atb[rr];
    }

    const double tol = 1e-10;
    for(int col = 0; col < dim_; col++){
        int pivot = col;
        for(int row = col + 1; row < dim_; row++){
            if(std::fabs(aug[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]) >
               std::fabs(aug[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(col)]))
                pivot = row;
        }
        if(std::fabs(aug[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(col)]) <= tol)
            return weighted_barycenter(anchors);
        if(pivot != col)
            std::swap(aug[static_cast<std::size_t>(col)], aug[static_cast<std::size_t>(pivot)]);

        const double pivot_val = aug[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)];
        for(int k = col; k <= dim_; k++)
            aug[static_cast<std::size_t>(col)][static_cast<std::size_t>(k)] /= pivot_val;
        for(int row = 0; row < dim_; row++){
            if(row == col)
                continue;
            const double factor = aug[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
            if(factor == 0.0)
                continue;
            for(int k = col; k <= dim_; k++){
                aug[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)] -=
                    factor * aug[static_cast<std::size_t>(col)][static_cast<std::size_t>(k)];
            }
        }
    }

    PointND out(static_cast<std::size_t>(dim_));
    for(int k = 0; k < dim_; k++)
        out[static_cast<std::size_t>(k)] =
            aug[static_cast<std::size_t>(k)][static_cast<std::size_t>(dim_)];
    return out;
}

PointND DrawGraphND::weighted_circle_position(
    const std::vector<std::pair<vertex_t, double>> &anchors) const
{
    if(dim_ != 2 || anchors.size() < 3)
        return weighted_barycenter(anchors);

    auto add_intersections = [&](const PointND &p1, double r1,
                                 const PointND &p2, double r2,
                                 std::vector<PointND> &out) {
        const double dx = p2[0] - p1[0];
        const double dy = p2[1] - p1[1];
        const double d = std::sqrt(dx * dx + dy * dy);
        if(d <= 0.0 || d > r1 + r2 || d < std::fabs(r1 - r2))
            return;
        const double a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
        double h2 = r1 * r1 - a * a;
        if(h2 < 0.0)
            h2 = 0.0;
        const double h = std::sqrt(h2);
        const double xm = p1[0] + a * dx / d;
        const double ym = p1[1] + a * dy / d;
        const double rx = -dy * (h / d);
        const double ry = dx * (h / d);

        PointND first(2);
        first[0] = xm + rx;
        first[1] = ym + ry;
        out.push_back(first);
        if(h2 > 0.0){
            PointND second(2);
            second[0] = xm - rx;
            second[1] = ym - ry;
            out.push_back(second);
        }
    };

    std::vector<PointND> candidates;
    const double r0 = anchors[0].second * kLegacyEdgeND;
    const double r1 = anchors[1].second * kLegacyEdgeND;
    const double r2 = anchors[2].second * kLegacyEdgeND;
    add_intersections(coords_[anchors[0].first], r0, coords_[anchors[1].first], r1, candidates);
    add_intersections(coords_[anchors[1].first], r1, coords_[anchors[2].first], r2, candidates);
    add_intersections(coords_[anchors[0].first], r0, coords_[anchors[2].first], r2, candidates);
    if(candidates.size() < 3)
        return weighted_barycenter(anchors);

    double best_score = std::numeric_limits<double>::infinity();
    PointND best(2);
    for(std::size_t i = 0; i < candidates.size(); i++){
        for(std::size_t j = i + 1; j < candidates.size(); j++){
            for(std::size_t k = j + 1; k < candidates.size(); k++){
                const double score =
                    std::sqrt(squared_distance(candidates[i], candidates[j])) +
                    std::sqrt(squared_distance(candidates[i], candidates[k])) +
                    std::sqrt(squared_distance(candidates[j], candidates[k]));
                if(score < best_score){
                    best_score = score;
                    best[0] = (candidates[i][0] + candidates[j][0] + candidates[k][0]) / 3.0;
                    best[1] = (candidates[i][1] + candidates[j][1] + candidates[k][1]) / 3.0;
                }
            }
        }
    }
    return std::isfinite(best_score) ? best : weighted_barycenter(anchors);
}

PointND DrawGraphND::local_kk_displacement(
    vertex_t root,
    const std::vector<std::pair<vertex_t, double>> &anchors) const
{
    PointND disp(static_cast<std::size_t>(dim_));
    for(const auto &anchor : anchors){
        if(anchor.second <= 0.0)
            continue;
        PointND vect(static_cast<std::size_t>(dim_));
        double norm2 = 0.0;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] = coords_[anchor.first][dd] - coords_[root][dd];
            norm2 += vect[dd] * vect[dd];
        }
        const double target2 = anchor.second * anchor.second * kLegacyEdgeND * kLegacyEdgeND;
        if(target2 <= 0.0)
            continue;
        const double scale = norm2 / target2 - 1.0;
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            vect[dd] *= scale;
            disp[dd] += vect[dd];
        }
    }
    return disp;
}

void DrawGraphND::polish_inserted_vertex(
    vertex_t root,
    const std::vector<std::pair<vertex_t, double>> &anchors,
    int steps)
{
    if(anchors.empty() || steps <= 0)
        return;

    for(int iter = 0; iter < steps; iter++){
        disp_[root] = local_kk_displacement(root, anchors);
        disp_norm_[root] = disp_[root].norm();
        if(disp_norm_[root] <= 0.0 || !std::isfinite(disp_norm_[root]))
            return;
        update_local_temperature(root);
        old_disp_[root] = disp_[root];
        old_disp_norm_[root] = disp_norm_[root];
        for(int d = 0; d < dim_; d++)
            disp_[root][static_cast<std::size_t>(d)] *= heat_[root];
        if(disp_norm_[root] > 0.0){
            for(int d = 0; d < dim_; d++)
                disp_[root][static_cast<std::size_t>(d)] /= disp_norm_[root];
        }
        for(int d = 0; d < dim_; d++){
            const std::size_t dd = static_cast<std::size_t>(d);
            coords_[root][dd] += disp_[root][dd];
            if(!std::isfinite(coords_[root][dd]))
                coords_[root][dd] = 0.0;
        }
    }
}

std::vector<std::pair<vertex_t, double>> DrawGraphND::weighted_anchor_distances(
    vertex_t root,
    const WeightedMisfND &misf,
    int active_count,
    int root_depth,
    int target_count,
    int strategy) const
{
    std::vector<char> eligible(graph_.size(), 0);
    for(int i = 0; i < active_count && i < static_cast<int>(misf.order.size()); i++){
        const vertex_t vert = misf.order[static_cast<std::size_t>(i)];
        if(insertion_anchor_scope_ == INSERT_ANCHOR_SCOPE_PREV_MISF_ND){
            if(misf.vertex_depth[vert] == root_depth + 1)
                eligible[vert] = 1;
        } else if(misf.vertex_depth[vert] > root_depth) {
            eligible[vert] = 1;
        }
    }

    const double inf = std::numeric_limits<double>::infinity();
    struct Node {
        double dist;
        vertex_t vert;
    };
    struct Greater {
        bool operator()(const Node &lhs, const Node &rhs) const
        {
            if(lhs.dist != rhs.dist)
                return lhs.dist > rhs.dist;
            return lhs.vert > rhs.vert;
        }
    };

    std::vector<double> dist(graph_.size(), inf);
    std::priority_queue<Node, std::vector<Node>, Greater> pq;
    dist[root] = 0.0;
    pq.push(Node{0.0, root});
    std::vector<std::pair<vertex_t, double>> anchors;
    double cutoff_dist = 0.0;
    int settled = 0;

    while(!pq.empty()){
        Node node = pq.top();
        pq.pop();
        if(node.dist > dist[node.vert] + kMetricTolND)
            continue;
        if(node.vert != root){
            settled++;
            if(eligible[node.vert]){
                anchors.push_back(std::make_pair(node.vert, node.dist));
                if(strategy == INSERT_ANCHOR_STRATEGY_FIRST_ND &&
                   static_cast<int>(anchors.size()) >= target_count)
                    break;
                if(strategy != INSERT_ANCHOR_STRATEGY_FIRST_ND &&
                   static_cast<int>(anchors.size()) == target_count)
                    cutoff_dist = node.dist;
            }
            if(strategy != INSERT_ANCHOR_STRATEGY_FIRST_ND &&
               cutoff_dist > 0.0 &&
               node.dist > cutoff_dist + kMetricTolND)
                break;
            if(metric_neighbor_cap_ > 0 && settled >= metric_neighbor_cap_)
                break;
        }
        for(const NeighborND &nbr : graph_.neighbors(node.vert)){
            const double alt = node.dist + nbr.weight;
            const double best = dist[nbr.vertex];
            const double scale = std::max(1.0,
                                          std::max(std::fabs(alt),
                                                   std::isfinite(best) ? std::fabs(best) : 0.0));
            if(!std::isfinite(best) || alt + kMetricTolND * scale < best){
                dist[nbr.vertex] = alt;
                pq.push(Node{alt, nbr.vertex});
            }
        }
    }
    if(strategy == INSERT_ANCHOR_STRATEGY_FIRST_ND)
        return anchors;
    return select_insertion_anchor_subset(anchors, target_count);
}

double DrawGraphND::weighted_eccentricity(vertex_t root) const
{
    const double inf = std::numeric_limits<double>::infinity();
    struct Node {
        double dist;
        vertex_t vert;
    };
    struct Greater {
        bool operator()(const Node &lhs, const Node &rhs) const
        {
            if(lhs.dist != rhs.dist)
                return lhs.dist > rhs.dist;
            return lhs.vert > rhs.vert;
        }
    };

    std::vector<double> dist(graph_.size(), inf);
    std::priority_queue<Node, std::vector<Node>, Greater> pq;
    dist[root] = 0.0;
    pq.push(Node{0.0, root});
    double out = 0.0;

    while(!pq.empty()){
        Node node = pq.top();
        pq.pop();
        if(node.dist > dist[node.vert] + kMetricTolND)
            continue;
        out = std::max(out, node.dist);
        for(const NeighborND &nbr : graph_.neighbors(node.vert)){
            const double alt = node.dist + nbr.weight;
            const double best = dist[nbr.vertex];
            const double scale = std::max(1.0,
                                          std::max(std::fabs(alt),
                                                   std::isfinite(best) ? std::fabs(best) : 0.0));
            if(!std::isfinite(best) || alt + kMetricTolND * scale < best){
                dist[nbr.vertex] = alt;
                pq.push(Node{alt, nbr.vertex});
            }
        }
    }
    return out;
}

int DrawGraphND::scheduled_rounds(int active_count) const
{
    if(active_count <= 0)
        return final_rounds_;
    if(active_count >= static_cast<int>(graph_.size()))
        return final_rounds_;
    if(rounds_ <= 0 || final_rounds_ <= 0)
        return std::max(0, final_rounds_);

    const double k = -std::log(static_cast<double>(final_rounds_) /
                               static_cast<double>(rounds_)) /
        static_cast<double>(std::max<std::size_t>(1, graph_.size()));
    return static_cast<int>(std::ceil(static_cast<double>(rounds_) *
                                      std::exp(-k * static_cast<double>(active_count))));
}

int DrawGraphND::metric_neighbor_cap_for_root(
    vertex_t root,
    const WeightedMisfND &misf) const
{
    if(metric_neighbor_cap_ <= 0 ||
       root >= misf.inverse.size() ||
       static_cast<int>(misf.inverse[root]) < misf.num_init)
        return 0;
    return metric_neighbor_cap_;
}

unsigned long DrawGraphND::next_legacy_rand()
{
    legacy_rng_state_ = 1664525UL * legacy_rng_state_ + 1013904223UL;
    return legacy_rng_state_;
}

PointND DrawGraphND::random_legacy_point()
{
    PointND out(static_cast<std::size_t>(dim_));
    const unsigned long box2 =
        static_cast<unsigned long>(std::max(1, box2_size_));
    for(int d = 0; d < dim_; d++)
        out[static_cast<std::size_t>(d)] =
            static_cast<double>(next_legacy_rand() % box2) - box_size_;
    return out;
}

bool DrawGraphND::should_record_refinement_step(int level_index,
                                                int misf_level,
                                                int round) const
{
    return refinement_step_trace_enabled_ &&
        level_index == refinement_step_trace_level_index_ &&
        misf_level == refinement_step_trace_misf_level_ &&
        round >= refinement_step_trace_round_start_ &&
        round <= refinement_step_trace_round_end_;
}

std::size_t DrawGraphND::record_refinement_step_pre(
    vertex_t vert,
    int order_index,
    int active_count,
    int level_index,
    int misf_level,
    int round,
    const PointND &coord_before,
    const PointND &pre_temp_disp,
    const std::vector<double> &attraction_disp,
    const std::vector<double> &repulsion_disp,
    const PointND &applied_disp,
    double heat_before,
    double heat_after,
    double old_cos_before,
    double old_cos_after,
    double old_disp_norm_before,
    double pre_temp_disp_norm,
    const std::string &attraction_edges,
    const std::string &repulsion_neighbors)
{
    std::vector<double> coord_before_flat(static_cast<std::size_t>(dim_));
    std::vector<double> pre_temp_disp_flat(static_cast<std::size_t>(dim_));
    std::vector<double> attraction_disp_flat(static_cast<std::size_t>(dim_));
    std::vector<double> repulsion_disp_flat(static_cast<std::size_t>(dim_));
    std::vector<double> applied_disp_flat(static_cast<std::size_t>(dim_));
    std::vector<double> coord_after_flat(static_cast<std::size_t>(dim_),
                                         std::numeric_limits<double>::quiet_NaN());
    for(int d = 0; d < dim_; d++){
        const std::size_t dd = static_cast<std::size_t>(d);
        coord_before_flat[dd] = coord_before[dd];
        pre_temp_disp_flat[dd] = pre_temp_disp[dd];
        attraction_disp_flat[dd] = attraction_disp[dd];
        repulsion_disp_flat[dd] = repulsion_disp[dd];
        applied_disp_flat[dd] = applied_disp[dd];
    }

    refinement_step_trace_.level_indices.push_back(level_index);
    refinement_step_trace_.misf_levels.push_back(misf_level);
    refinement_step_trace_.rounds.push_back(round);
    refinement_step_trace_.active_counts.push_back(active_count);
    refinement_step_trace_.order_indices.push_back(order_index);
    refinement_step_trace_.vertices.push_back(static_cast<int>(vert) + 1);
    refinement_step_trace_.heat_before.push_back(heat_before);
    refinement_step_trace_.heat_after.push_back(heat_after);
    refinement_step_trace_.old_cos_before.push_back(old_cos_before);
    refinement_step_trace_.old_cos_after.push_back(old_cos_after);
    refinement_step_trace_.old_disp_norm_before.push_back(old_disp_norm_before);
    refinement_step_trace_.pre_temp_disp_norm.push_back(pre_temp_disp_norm);
    refinement_step_trace_.coords_before.push_back(std::move(coord_before_flat));
    refinement_step_trace_.pre_temp_disp.push_back(std::move(pre_temp_disp_flat));
    refinement_step_trace_.attraction_disp.push_back(std::move(attraction_disp_flat));
    refinement_step_trace_.repulsion_disp.push_back(std::move(repulsion_disp_flat));
    refinement_step_trace_.applied_disp.push_back(std::move(applied_disp_flat));
    refinement_step_trace_.coords_after.push_back(std::move(coord_after_flat));
    refinement_step_trace_.attraction_edges.push_back(attraction_edges);
    refinement_step_trace_.repulsion_neighbors.push_back(repulsion_neighbors);
    const int parent_row = static_cast<int>(refinement_step_trace_.vertices.size());
    for(std::size_t i = 0; i < last_attraction_term_neighbors_.size(); i++){
        refinement_step_trace_.attraction_term_parent_rows.push_back(parent_row);
        refinement_step_trace_.attraction_term_indices.push_back(static_cast<int>(i) + 1);
        refinement_step_trace_.attraction_term_vertices.push_back(static_cast<int>(vert) + 1);
        refinement_step_trace_.attraction_term_neighbors.push_back(
            last_attraction_term_neighbors_[i]);
        refinement_step_trace_.attraction_term_weights.push_back(last_attraction_term_weights_[i]);
        refinement_step_trace_.attraction_term_norm2.push_back(last_attraction_term_norm2_[i]);
        refinement_step_trace_.attraction_term_desired.push_back(last_attraction_term_desired_[i]);
        refinement_step_trace_.attraction_term_desired2.push_back(
            last_attraction_term_desired2_[i]);
        refinement_step_trace_.attraction_term_scale.push_back(last_attraction_term_scale_[i]);
        refinement_step_trace_.attraction_term_delta.push_back(last_attraction_term_delta_[i]);
        refinement_step_trace_.attraction_term_step.push_back(last_attraction_term_step_[i]);
        refinement_step_trace_.attraction_term_cumulative.push_back(
            last_attraction_term_cumulative_[i]);
    }
    return refinement_step_trace_.vertices.size() - 1;
}

void DrawGraphND::record_refinement_step_after(std::size_t row, vertex_t vert)
{
    if(row >= refinement_step_trace_.coords_after.size())
        return;
    for(int d = 0; d < dim_; d++)
        refinement_step_trace_.coords_after[row][static_cast<std::size_t>(d)] =
            coords_[vert][static_cast<std::size_t>(d)];
}

void DrawGraphND::record_insertion_trace(
    vertex_t root,
    int level_index,
    int misf_level,
    int previous_active_count,
    int active_count,
    int order_index,
    int root_depth,
    int anchor_count_requested,
    int anchor_count_used,
    int insertion_mode,
    int local_kk_steps,
    const std::vector<std::pair<vertex_t, double>> &anchors,
    const PointND &coord_initial,
    const PointND &coord_after,
    const PointND &old_disp_initial,
    const PointND &old_disp_after,
    double old_disp_norm_initial_arg,
    double old_disp_norm_after_arg)
{
    if(!insertion_trace_enabled_)
        return;

    auto flatten_point = [&](const PointND &point) {
        std::vector<double> out(static_cast<std::size_t>(dim_), 0.0);
        for(int d = 0; d < dim_; d++)
            out[static_cast<std::size_t>(d)] = point[static_cast<std::size_t>(d)];
        return out;
    };

    std::ostringstream anchor_stream;
    anchor_stream << std::setprecision(17);
    const int limit = std::min(anchor_count_used, static_cast<int>(anchors.size()));
    for(int i = 0; i < limit; i++){
        if(i > 0)
            anchor_stream << ";";
        anchor_stream << static_cast<int>(anchors[static_cast<std::size_t>(i)].first) + 1
                      << ":"
                      << anchors[static_cast<std::size_t>(i)].second;
    }

    insertion_trace_.level_indices.push_back(level_index);
    insertion_trace_.misf_levels.push_back(misf_level);
    insertion_trace_.previous_active_counts.push_back(previous_active_count);
    insertion_trace_.active_counts.push_back(active_count);
    insertion_trace_.order_indices.push_back(order_index);
    insertion_trace_.vertices.push_back(static_cast<int>(root) + 1);
    insertion_trace_.root_depths.push_back(root_depth);
    insertion_trace_.anchor_count_requested.push_back(anchor_count_requested);
    insertion_trace_.anchor_count_used.push_back(anchor_count_used);
    insertion_trace_.insertion_modes.push_back(insertion_mode);
    insertion_trace_.local_kk_steps.push_back(local_kk_steps);
    insertion_trace_.anchors.push_back(anchor_stream.str());
    insertion_trace_.old_disp_norm_initial.push_back(old_disp_norm_initial_arg);
    insertion_trace_.old_disp_norm_after.push_back(old_disp_norm_after_arg);
    insertion_trace_.coords_initial.push_back(flatten_point(coord_initial));
    insertion_trace_.coords_after.push_back(flatten_point(coord_after));
    insertion_trace_.old_disp_initial.push_back(flatten_point(old_disp_initial));
    insertion_trace_.old_disp_after.push_back(flatten_point(old_disp_after));
}

} // namespace gripnd
