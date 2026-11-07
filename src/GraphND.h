#ifndef GRIP_GRAPH_ND_H
#define GRIP_GRAPH_ND_H

#include <cstdint>
#include <vector>

namespace gripnd {

typedef std::uint32_t vertex_t;

struct NeighborND {
    vertex_t vertex;
    double weight;
};

struct EdgeND {
    vertex_t source;
    vertex_t target;
    double weight;
};

class GraphND {
public:
    explicit GraphND(std::size_t n = 0);

    std::size_t size() const { return adj_.size(); }
    void add_directed_edge(vertex_t source, vertex_t target, double weight);

    const std::vector<NeighborND> &neighbors(vertex_t vertex) const
    {
        return adj_[vertex];
    }

    std::vector<EdgeND> canonical_edges() const;
    double median_edge_weight() const;

private:
    std::vector<std::vector<NeighborND>> adj_;
};

} // namespace gripnd

#endif
