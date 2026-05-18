#ifndef GRIP_WEIGHTED_MISF_ND_H
#define GRIP_WEIGHTED_MISF_ND_H

#include "GraphND.h"

#include <vector>

namespace gripnd {

struct WeightedMisfND {
    std::vector<vertex_t> order;
    std::vector<vertex_t> inverse;
    std::vector<int> vertex_depth;
    std::vector<int> level_size;
    std::vector<int> num_nbrs_schedule;
    int height;
    int num_init;
    unsigned long rng_state;
};

WeightedMisfND build_weighted_misf_nd(const GraphND &graph,
                                      int num_init,
                                      int num_nbrs,
                                      unsigned long seed);

} // namespace gripnd

#endif
