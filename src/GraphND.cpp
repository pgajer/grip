#include "GraphND.h"

#include <algorithm>
#include <set>

namespace gripnd {

GraphND::GraphND(std::size_t n) : adj_(n) {}

void GraphND::add_directed_edge(vertex_t source, vertex_t target, double weight)
{
    adj_[source].push_back(NeighborND{target, weight});
}

std::vector<EdgeND> GraphND::canonical_edges() const
{
    std::vector<EdgeND> edges;
    std::set<std::pair<vertex_t, vertex_t>> seen;
    for(vertex_t source = 0; source < adj_.size(); source++){
        for(const NeighborND &nbr : adj_[source]){
            const vertex_t a = std::min(source, nbr.vertex);
            const vertex_t b = std::max(source, nbr.vertex);
            if(a == b)
                continue;
            if(seen.insert(std::make_pair(a, b)).second)
                edges.push_back(EdgeND{a, b, nbr.weight});
        }
    }
    return edges;
}

double GraphND::median_edge_weight() const
{
    std::vector<double> weights;
    std::set<std::pair<vertex_t, vertex_t>> seen;
    for(vertex_t source = 0; source < adj_.size(); source++){
        for(const NeighborND &nbr : adj_[source]){
            const vertex_t a = std::min(source, nbr.vertex);
            const vertex_t b = std::max(source, nbr.vertex);
            if(a != b && seen.insert(std::make_pair(a, b)).second)
                weights.push_back(nbr.weight);
        }
    }
    if(weights.empty())
        return 1.0;

    const std::size_t mid = weights.size() / 2;
    std::nth_element(weights.begin(), weights.begin() + mid, weights.end());
    double med = weights[mid];
    if(weights.size() % 2 == 0){
        std::nth_element(weights.begin(), weights.begin() + mid - 1, weights.end());
        med = 0.5 * (med + weights[mid - 1]);
    }
    return med > 0.0 ? med : 1.0;
}

} // namespace gripnd
