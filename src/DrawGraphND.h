#ifndef GRIP_DRAW_GRAPH_ND_H
#define GRIP_DRAW_GRAPH_ND_H

#include "GraphND.h"
#include "PointND.h"

#include <random>
#include <vector>

namespace gripnd {

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
                unsigned int seed);

    std::vector<PointND> layout();

private:
    void initialize();
    void refine(const std::vector<EdgeND> &edges, int total_rounds);
    void add_edge_forces(const std::vector<EdgeND> &edges,
                         std::vector<PointND> &disp,
                         double edge_rate) const;
    void add_repulsion_forces(std::vector<PointND> &disp, double repulse_rate) const;
    void apply_displacements(const std::vector<PointND> &disp, double max_move);

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
    std::mt19937 rng_;
    std::vector<PointND> coords_;
};

} // namespace gripnd

#endif
