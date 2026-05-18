#include "DrawGraphND.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace gripnd {

DrawGraphND::DrawGraphND(const GraphND &graph,
                         int dim,
                         int rounds,
                         int final_rounds,
                         int num_init,
                         int num_nbrs,
                         double r,
                         double s,
                         double repulsion_factor,
                         int tinit_factor,
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
      tinit_factor_(tinit_factor),
      rng_(seed),
      coords_(graph.size(), PointND(static_cast<std::size_t>(dim)))
{
}

std::vector<PointND> DrawGraphND::layout()
{
    initialize();
    refine(graph_.canonical_edges(), rounds_ + final_rounds_);
    return coords_;
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

} // namespace gripnd
