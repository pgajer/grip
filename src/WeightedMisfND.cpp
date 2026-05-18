#include "WeightedMisfND.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>

namespace gripnd {
namespace {

constexpr double kMetricTolND = 1e-10;

int ilog2_floor(std::size_t n)
{
    int k = 0;
    while(n > 0){
        n /= 2;
        k++;
    }
    return k - 1;
}

double weighted_misf_radius_nd(int level)
{
    if(level <= 1)
        return 1.0;
    return std::pow(2.0, static_cast<double>(level - 1));
}

class FastRandND {
public:
    explicit FastRandND(unsigned long seed) : state_(seed) {}

    unsigned long next()
    {
        state_ = 1664525UL * state_ + 1013904223UL;
        return state_;
    }

    unsigned long state() const { return state_; }

private:
    unsigned long state_;
};

template <class Visitor>
void traverse_weighted_shortest_paths_nd(const GraphND &graph,
                                         vertex_t root,
                                         double cutoff,
                                         Visitor visitor)
{
    const double inf = std::numeric_limits<double>::infinity();
    if(root >= graph.size())
        return;

    struct QueueNode {
        double dist;
        vertex_t vert;
    };
    struct QueueNodeGreater {
        bool operator()(const QueueNode &lhs, const QueueNode &rhs) const
        {
            if(lhs.dist != rhs.dist)
                return lhs.dist > rhs.dist;
            return lhs.vert > rhs.vert;
        }
    };

    std::vector<double> dist(graph.size(), inf);
    std::priority_queue<QueueNode,
                        std::vector<QueueNode>,
                        QueueNodeGreater> pq;
    dist[root] = 0.0;
    pq.push(QueueNode{0.0, root});

    while(!pq.empty()){
        QueueNode node = pq.top();
        pq.pop();
        if(node.dist > dist[node.vert] + kMetricTolND)
            continue;
        if(node.dist > cutoff + kMetricTolND)
            break;

        if(node.vert != root){
            if(visitor(node.vert, node.dist))
                break;
        }

        const std::vector<NeighborND> &nbrs = graph.neighbors(node.vert);
        for(const NeighborND &nbr : nbrs){
            double alt = node.dist + nbr.weight;
            if(alt > cutoff + kMetricTolND)
                continue;
            double best = dist[nbr.vertex];
            double scale = std::max(1.0,
                                    std::max(std::fabs(alt),
                                             std::isfinite(best) ? std::fabs(best) : 0.0));
            if(!std::isfinite(best) || alt + kMetricTolND * scale < best){
                dist[nbr.vertex] = alt;
                pq.push(QueueNode{alt, nbr.vertex});
            }
        }
    }
}

} // namespace

WeightedMisfND build_weighted_misf_nd(const GraphND &graph,
                                      int num_init,
                                      int num_nbrs,
                                      unsigned long seed)
{
    const int n = static_cast<int>(graph.size());
    num_init = std::max(1, std::min(num_init, n));
    num_nbrs = std::max(1, num_nbrs);

    WeightedMisfND out;
    out.order.resize(graph.size());
    out.inverse.resize(graph.size());
    out.vertex_depth.assign(graph.size(), 0);
    out.height = 0;
    out.num_init = num_init;
    out.rng_state = seed;

    for(int i = 0; i < n; i++){
        out.order[static_cast<std::size_t>(i)] = static_cast<vertex_t>(i);
        out.inverse[static_cast<std::size_t>(i)] = static_cast<vertex_t>(i);
    }

    const int log_2_n = ilog2_floor(graph.size()) + 2;
    std::vector<int> misf_size(static_cast<std::size_t>(log_2_n), 0);

    if(n == num_init){
        misf_size[0] = n;
        out.height = 0;
    } else {
        FastRandND rng(seed);
        int misf_level = 1;
        misf_size[0] = n;
        int curr_size = n;
        int itr = 0;

        do {
            const int prev_size = curr_size;
            curr_size = 0;
            itr = 0;
            const double radius = weighted_misf_radius_nd(misf_level);

            while(prev_size > num_init && itr < prev_size){
                const int vert =
                    itr + static_cast<int>(rng.next() % static_cast<unsigned long>(prev_size - itr));
                std::swap(out.order[static_cast<std::size_t>(vert)],
                          out.order[static_cast<std::size_t>(itr)]);
                out.inverse[out.order[static_cast<std::size_t>(vert)]] =
                    static_cast<vertex_t>(vert);
                out.inverse[out.order[static_cast<std::size_t>(itr)]] =
                    static_cast<vertex_t>(itr);

                std::swap(out.order[static_cast<std::size_t>(curr_size)],
                          out.order[static_cast<std::size_t>(itr)]);
                out.inverse[out.order[static_cast<std::size_t>(curr_size)]] =
                    static_cast<vertex_t>(curr_size);
                out.inverse[out.order[static_cast<std::size_t>(itr)]] =
                    static_cast<vertex_t>(itr);
                itr++;

                const vertex_t new_el = out.order[static_cast<std::size_t>(curr_size++)];
                out.vertex_depth[new_el] = misf_level;

                traverse_weighted_shortest_paths_nd(
                    graph,
                    new_el,
                    radius,
                    [&](vertex_t adj, double dist_adj) {
                        if(dist_adj > radius + kMetricTolND)
                            return true;
                        const int adj_inv = static_cast<int>(out.inverse[adj]);
                        if(adj_inv >= itr && adj_inv < prev_size){
                            std::swap(out.order[static_cast<std::size_t>(itr)],
                                      out.order[static_cast<std::size_t>(adj_inv)]);
                            out.inverse[out.order[static_cast<std::size_t>(adj_inv)]] =
                                static_cast<vertex_t>(adj_inv);
                            out.inverse[out.order[static_cast<std::size_t>(itr)]] =
                                static_cast<vertex_t>(itr);
                            itr++;
                        }
                        return false;
                    });
            }

            if(misf_level >= static_cast<int>(misf_size.size()))
                misf_size.resize(misf_size.size() + 1, 0);
            misf_size[static_cast<std::size_t>(misf_level)] = curr_size;
            misf_level++;
        } while(itr);

        misf_level -= 1;
        while(misf_level > 0 && misf_size[static_cast<std::size_t>(misf_level)] < num_init)
            misf_level--;

        if(misf_size[static_cast<std::size_t>(misf_level)] > num_init){
            int v = 0;
            while(v < num_init){
                out.vertex_depth[out.order[static_cast<std::size_t>(v)]] = misf_level + 1;
                v++;
            }
            misf_level++;
            if(misf_level >= static_cast<int>(misf_size.size()))
                misf_size.resize(misf_size.size() + 1, 0);
            misf_size[static_cast<std::size_t>(misf_level)] = num_init;
        }

        out.height = misf_level;
        out.rng_state = rng.state();
    }

    out.level_size.assign(misf_size.begin(), misf_size.begin() + out.height + 1);
    out.num_nbrs_schedule.resize(out.level_size.size());
    for(std::size_t level = 0; level < out.level_size.size(); level++){
        const int available = std::max(0, out.level_size[level] - 1);
        out.num_nbrs_schedule[level] = std::min(num_nbrs, available);
    }

    return out;
}

} // namespace gripnd
