#ifndef GRIP_DRAW_GRAPH_ND_H
#define GRIP_DRAW_GRAPH_ND_H

#include "GraphND.h"
#include "PointND.h"
#include "WeightedMisfND.h"

#include <random>
#include <string>
#include <vector>

namespace gripnd {

enum PlacementModeND {
    PLACEMENT_BARYCENTER_ND = 0,
    PLACEMENT_CIRCLE_ND = 1
};

enum InsertionAnchorScopeND {
    INSERT_ANCHOR_SCOPE_ANY_HIGHER_ND = 0,
    INSERT_ANCHOR_SCOPE_PREV_MISF_ND = 1
};

enum InsertionAnchorStrategyND {
    INSERT_ANCHOR_STRATEGY_FIRST_ND = 0,
    INSERT_ANCHOR_STRATEGY_DISTANCE_BAND_ND = 1,
    INSERT_ANCHOR_STRATEGY_BALANCED_BAND_ND = 2,
    INSERT_ANCHOR_STRATEGY_SPREAD_PREV_ND = 3
};

enum Level0InsertionModeND {
    LEVEL0_INSERT_INHERIT_ND = 100,
    LEVEL0_INSERT_BARYCENTER_ND = 101,
    LEVEL0_INSERT_LEAST_SQUARES_ND = 102
};

enum FinalStageModeND {
    FINAL_STAGE_FR_ND = 0,
    FINAL_STAGE_KK_REPULSE_ND = 1
};

struct LayoutTraceND {
    std::vector<std::vector<double>> frames;
    std::vector<std::string> phases;
    std::vector<int> rounds;
    std::vector<int> level_indices;
    std::vector<int> misf_levels;
    std::vector<int> active_counts;
};

struct RefinementStepTraceND {
    std::vector<int> level_indices;
    std::vector<int> misf_levels;
    std::vector<int> rounds;
    std::vector<int> active_counts;
    std::vector<int> order_indices;
    std::vector<int> vertices;
    std::vector<double> heat_before;
    std::vector<double> heat_after;
    std::vector<double> old_cos_before;
    std::vector<double> old_cos_after;
    std::vector<double> old_disp_norm_before;
    std::vector<double> pre_temp_disp_norm;
    std::vector<std::vector<double>> coords_before;
    std::vector<std::vector<double>> pre_temp_disp;
    std::vector<std::vector<double>> attraction_disp;
    std::vector<std::vector<double>> repulsion_disp;
    std::vector<std::vector<double>> applied_disp;
    std::vector<std::vector<double>> coords_after;
    std::vector<std::string> attraction_edges;
    std::vector<std::string> repulsion_neighbors;
    std::vector<int> attraction_term_parent_rows;
    std::vector<int> attraction_term_indices;
    std::vector<int> attraction_term_vertices;
    std::vector<int> attraction_term_neighbors;
    std::vector<double> attraction_term_weights;
    std::vector<double> attraction_term_norm2;
    std::vector<double> attraction_term_desired;
    std::vector<double> attraction_term_desired2;
    std::vector<double> attraction_term_scale;
    std::vector<std::vector<double>> attraction_term_delta;
    std::vector<std::vector<double>> attraction_term_step;
    std::vector<std::vector<double>> attraction_term_cumulative;
};

struct InsertionTraceND {
    std::vector<int> level_indices;
    std::vector<int> misf_levels;
    std::vector<int> previous_active_counts;
    std::vector<int> active_counts;
    std::vector<int> order_indices;
    std::vector<int> vertices;
    std::vector<int> root_depths;
    std::vector<int> anchor_count_requested;
    std::vector<int> anchor_count_used;
    std::vector<int> insertion_modes;
    std::vector<int> local_kk_steps;
    std::vector<std::string> anchors;
    std::vector<double> old_disp_norm_initial;
    std::vector<double> old_disp_norm_after;
    std::vector<std::vector<double>> coords_initial;
    std::vector<std::vector<double>> coords_after;
    std::vector<std::vector<double>> old_disp_initial;
    std::vector<std::vector<double>> old_disp_after;
};

struct MetricNeighborND {
    vertex_t vert;
    double dist;
};

class DrawGraphND {
public:
    DrawGraphND(const GraphND &graph,
                int dim,
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
                int final_stage_mode,
                int metric_neighbor_cap,
                int placement_mode,
                int insertion_anchor_count,
                int insertion_anchor_scope,
                int insertion_anchor_strategy,
                int level0_insertion_mode,
                int level0_anchor_count,
                int level0_local_kk_steps,
                unsigned int seed);

    std::vector<PointND> layout();
    std::vector<PointND> layout(LayoutTraceND *trace, int trace_every);
    void configure_refinement_step_trace(int level_index,
                                         int misf_level,
                                         int round_start,
                                         int round_end);
    void configure_insertion_trace(bool enabled);
    const RefinementStepTraceND &get_refinement_step_trace() const
    {
        return refinement_step_trace_;
    }
    const InsertionTraceND &get_insertion_trace() const
    {
        return insertion_trace_;
    }

private:
    void initialize();
    void refine(const std::vector<EdgeND> &edges, int total_rounds);
    void refine_trace_level(const std::vector<EdgeND> &edges,
                            int level_rounds,
                            int active_count,
                            LayoutTraceND *trace,
                            int trace_every,
                            int level_index,
                            int misf_level);
    void add_edge_forces(const std::vector<EdgeND> &edges,
                         std::vector<PointND> &disp,
                         double edge_rate) const;
    void add_repulsion_forces(std::vector<PointND> &disp, double repulse_rate) const;
    void apply_displacements(const std::vector<PointND> &disp, double max_move);
    void apply_displacements(const std::vector<PointND> &disp,
                             double max_move,
                             const std::vector<char> &active);
    void record_trace(LayoutTraceND *trace,
                      const std::string &phase,
                      int round,
                      int active_count,
                      int level_index,
                      int misf_level) const;
    void reset_active_heat(int active_count);
    void refine_legacy_weighted_level(const WeightedMisfND &misf,
                                      int level_rounds,
                                      int active_count,
                                      LayoutTraceND *trace,
                                      int trace_every,
                                      int level_index,
                                      int misf_level);
    std::vector<MetricNeighborND> metric_neighbors_for_level(vertex_t root,
                                                             const WeightedMisfND &misf,
                                                             int target_level) const;
    void populate_metric_neighbors_for_root(vertex_t root,
                                            const WeightedMisfND &misf) const;
    void legacy_weighted_kk_displacement(vertex_t vert,
                                         const std::vector<MetricNeighborND> &neighbors,
                                         int active_count,
                                         int misf_level);
    void legacy_weighted_kk_final_displacement(vertex_t vert,
                                               const std::vector<MetricNeighborND> &neighbors,
                                               int active_count);
    void legacy_weighted_fr_displacement(vertex_t vert,
                                         const std::vector<MetricNeighborND> &neighbors,
                                         int active_count);
    void add_legacy_active_repulsion(vertex_t vert,
                                     int active_count,
                                     double repulsion_scale);
    void add_final_anchor_force(vertex_t vert);
    void prepare_final_anchors(int active_count);
    void scale_legacy_displacement(vertex_t vert);
    void update_local_temperature(vertex_t vert);
    void initialize_multiscale_trace(const WeightedMisfND &misf);
    void initialize_top_level(const WeightedMisfND &misf);
    void insert_level_vertices(const WeightedMisfND &misf,
                               int previous_active_count,
                               int active_count,
                               int level_index,
                               int misf_level);
    std::vector<std::pair<vertex_t, double>> select_insertion_anchor_subset(
        const std::vector<std::pair<vertex_t, double>> &anchors,
        int target_count) const;
    PointND weighted_barycenter(const std::vector<vertex_t> &anchors) const;
    PointND weighted_barycenter(const std::vector<std::pair<vertex_t, double>> &anchors) const;
    PointND weighted_least_squares_position(const std::vector<std::pair<vertex_t, double>> &anchors) const;
    PointND weighted_circle_position(const std::vector<std::pair<vertex_t, double>> &anchors) const;
    PointND local_kk_displacement(vertex_t root,
                                  const std::vector<std::pair<vertex_t, double>> &anchors) const;
    void polish_inserted_vertex(vertex_t root,
                                const std::vector<std::pair<vertex_t, double>> &anchors,
                                int steps);
    std::vector<std::pair<vertex_t, double>> weighted_anchor_distances(vertex_t root,
                                                                       const WeightedMisfND &misf,
                                                                       int active_count,
                                                                       int root_depth,
                                                                       int target_count,
                                                                       int strategy) const;
    double weighted_eccentricity(vertex_t root) const;
    int scheduled_rounds(int active_count) const;
    int metric_neighbor_cap_for_root(vertex_t root,
                                     const WeightedMisfND &misf) const;
    unsigned long next_legacy_rand();
    PointND random_legacy_point();
    bool should_record_refinement_step(int level_index,
                                       int misf_level,
                                       int round) const;
    std::size_t record_refinement_step_pre(vertex_t vert,
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
                                           const std::string &repulsion_neighbors);
    void record_refinement_step_after(std::size_t row, vertex_t vert);
    void record_insertion_trace(vertex_t root,
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
                                double old_disp_norm_initial,
                                double old_disp_norm_after);

    const GraphND &graph_;
    int dim_;
    int rounds_;
    int final_rounds_;
    int num_init_;
    int num_nbrs_;
    double r_;
    double s_;
    double repulsion_factor_;
    int tinit_factor_;
    double final_anchor_factor_;
    double final_move_scale_after_first_;
    int final_stage_mode_;
    int metric_neighbor_cap_;
    int placement_mode_;
    int insertion_anchor_count_;
    int insertion_anchor_scope_;
    int insertion_anchor_strategy_;
    int level0_insertion_mode_;
    int level0_anchor_count_;
    int level0_local_kk_steps_;
    unsigned int seed_;
    std::mt19937 rng_;
    unsigned long legacy_rng_state_;
    double box_size_;
    int box2_size_;
    std::vector<PointND> coords_;
    std::vector<PointND> disp_;
    std::vector<PointND> old_disp_;
    std::vector<double> disp_norm_;
    std::vector<double> old_disp_norm_;
    std::vector<double> heat_;
    std::vector<double> old_cos_;
    bool final_anchor_ready_;
    std::vector<PointND> final_anchor_pos_;
    std::vector<vertex_t> trace_order_;
    mutable std::vector<std::vector<std::vector<MetricNeighborND>>> metric_neighbors_cache_;
    mutable std::vector<unsigned char> metric_neighbors_cached_;
    bool refinement_step_trace_enabled_;
    int refinement_step_trace_level_index_;
    int refinement_step_trace_misf_level_;
    int refinement_step_trace_round_start_;
    int refinement_step_trace_round_end_;
    RefinementStepTraceND refinement_step_trace_;
    bool insertion_trace_enabled_;
    InsertionTraceND insertion_trace_;
    std::vector<double> last_attraction_disp_;
    std::vector<double> last_repulsion_disp_;
    std::string last_attraction_edges_;
    std::string last_repulsion_neighbors_;
    std::vector<int> last_attraction_term_neighbors_;
    std::vector<double> last_attraction_term_weights_;
    std::vector<double> last_attraction_term_norm2_;
    std::vector<double> last_attraction_term_desired_;
    std::vector<double> last_attraction_term_desired2_;
    std::vector<double> last_attraction_term_scale_;
    std::vector<std::vector<double>> last_attraction_term_delta_;
    std::vector<std::vector<double>> last_attraction_term_step_;
    std::vector<std::vector<double>> last_attraction_term_cumulative_;
};

} // namespace gripnd

#endif
