// Weighted GRIP support routines

#include "DrawGraph.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <vector>

namespace {

constexpr double kMetricTol = 1e-10;

double weighted_misf_radius(size_tt misfLevel)
{
    if(misfLevel <= 1)
        return 1.0;
    return std::pow(2.0, static_cast<double>(misfLevel - 1));
}

} // namespace

void DrawGraph::traverse_weighted_shortest_paths(
    size_tt root,
    double cutoff,
    size_tt maxSettled,
    const std::function<bool(size_tt, double)> &visitor)
{
    const double inf = std::numeric_limits<double>::infinity();
    if(root >= numOfVert)
        return;

    if(metricScratchDist.size() != numOfVert){
        metricScratchDist.resize(numOfVert, inf);
        metricScratchStamp.assign(numOfVert, 0);
        metricScratchEpoch = 0;
    }
    metricScratchEpoch++;
    if(metricScratchEpoch == 0){
        std::fill(metricScratchStamp.begin(), metricScratchStamp.end(), 0);
        metricScratchEpoch = 1;
    }

    struct QueueNode {
        double dist;
        size_tt vert;
    };
    struct QueueNodeGreater {
        bool operator()(const QueueNode &lhs, const QueueNode &rhs) const
        {
            if(lhs.dist != rhs.dist)
                return lhs.dist > rhs.dist;
            return lhs.vert > rhs.vert;
        }
    };

    auto get_dist = [&](size_tt vert){
        return metricScratchStamp[vert] == metricScratchEpoch
            ? metricScratchDist[vert]
            : inf;
    };
    auto set_dist = [&](size_tt vert, double value){
        metricScratchStamp[vert] = metricScratchEpoch;
        metricScratchDist[vert] = value;
    };

    std::priority_queue<QueueNode,
                        std::vector<QueueNode>,
                        QueueNodeGreater> pq;
    set_dist(root, 0.0);
    pq.push(QueueNode{0.0, root});

    size_tt settled = 0;
    while(!pq.empty()){
        QueueNode node = pq.top();
        pq.pop();

        double best = get_dist(node.vert);
        if(node.dist > best + kMetricTol)
            continue;
        if(node.dist > cutoff + kMetricTol)
            break;

        if(node.vert != root){
            settled++;
            if(visitor(node.vert, node.dist))
                break;
            if(maxSettled > 0 && settled >= maxSettled)
                break;
        }

        size_tt degLocal = graph.adjList[0][node.vert];
        for(size_tt adjVert = 0; adjVert < degLocal; adjVert++){
            size_tt overt = graph.adjList[node.vert + 1][adjVert];
            double alt = node.dist + graph.get_edge_weight(node.vert, adjVert);
            if(alt > cutoff + kMetricTol)
                continue;
            double overtBest = get_dist(overt);
            double scale = std::max(1.0,
                                    std::max(std::fabs(alt),
                                             std::isfinite(overtBest)
                                                 ? std::fabs(overtBest)
                                                 : 0.0));
            if(!std::isfinite(overtBest) || alt + kMetricTol * scale < overtBest){
                set_dist(overt, alt);
                pq.push(QueueNode{alt, overt});
            }
        }
    }
}

void DrawGraph::compute_weighted_shortest_paths(size_tt root,
                                                std::vector<double> &dist,
                                                double cutoff) const
{
    const double inf = std::numeric_limits<double>::infinity();
    dist.assign(numOfVert, inf);
    if(root >= numOfVert)
        return;

    struct QueueNode {
        double dist;
        size_tt vert;
    };
    struct QueueNodeGreater {
        bool operator()(const QueueNode &lhs, const QueueNode &rhs) const
        {
            if(lhs.dist != rhs.dist)
                return lhs.dist > rhs.dist;
            return lhs.vert > rhs.vert;
        }
    };

    std::priority_queue<QueueNode,
                        std::vector<QueueNode>,
                        QueueNodeGreater> pq;
    dist[root] = 0.0;
    pq.push(QueueNode{0.0, root});

    while(!pq.empty()){
        QueueNode node = pq.top();
        pq.pop();
        if(node.dist > dist[node.vert] + kMetricTol)
            continue;
        if(node.dist > cutoff + kMetricTol)
            continue;

        size_tt degLocal = graph.adjList[0][node.vert];
        for(size_tt adjVert = 0; adjVert < degLocal; adjVert++){
            size_tt overt = graph.adjList[node.vert + 1][adjVert];
            double alt = node.dist + graph.get_edge_weight(node.vert, adjVert);
            if(alt > cutoff + kMetricTol)
                continue;
            double best = dist[overt];
            double scale = std::max(1.0,
                                    std::max(std::fabs(alt),
                                             std::isfinite(best) ? std::fabs(best) : 0.0));
            if(!std::isfinite(best) || alt + kMetricTol * scale < best){
                dist[overt] = alt;
                pq.push(QueueNode{alt, overt});
            }
        }
    }
}

std::vector<size_tt> DrawGraph::ordered_vertices_by_metric(const std::vector<double> &dist,
                                                           size_tt root) const
{
    std::vector<size_tt> order;
    order.reserve(numOfVert > 0 ? numOfVert - 1 : 0);
    for(size_tt vert = 0; vert < numOfVert; vert++){
        if(vert == root || !std::isfinite(dist[vert]))
            continue;
        order.push_back(vert);
    }
    std::sort(order.begin(), order.end(),
              [&](size_tt lhs, size_tt rhs){
                  if(dist[lhs] != dist[rhs])
                      return dist[lhs] < dist[rhs];
                  return lhs < rhs;
              });
    return order;
}

void DrawGraph::create_misf_weighted()
{
    misfLevel = 1;
    misfSize[misfLevel - 1] = numOfVert;
    size_tt mishSizeCurrLevel = numOfVert;

    size_tt itr;
    do {
        size_tt mishSizePrevLevel = mishSizeCurrLevel;
        mishSizeCurrLevel = 0;
        itr = 0;
        double radius = weighted_misf_radius(misfLevel);

        while(mishSizePrevLevel > numOfInitVert && itr < mishSizePrevLevel){
            size_tt vert = itr + graph.fast_Rand() % (mishSizePrevLevel - itr);
            std::swap(mish[vert], mish[itr]);
            inv[mish[vert]] = vert;
            inv[mish[itr]] = itr;

            std::swap(mish[mishSizeCurrLevel], mish[itr]);
            inv[mish[mishSizeCurrLevel]] = mishSizeCurrLevel;
            inv[mish[itr]] = itr;
            itr++;

            size_tt newEl = mish[mishSizeCurrLevel++];
            vertDepth[newEl] = misfLevel;

            traverse_weighted_shortest_paths(
                newEl,
                radius,
                0,
                [&](size_tt adj, double distAdj){
                    if(distAdj > radius + kMetricTol)
                        return true;
                    if(inv[adj] >= itr && inv[adj] < mishSizePrevLevel){
                        std::swap(mish[itr], mish[inv[adj]]);
                        inv[mish[inv[adj]]] = inv[adj];
                        inv[mish[itr]] = itr;
                        itr++;
                    }
                    return false;
                });
        }

        if(misfLevel >= log_2_n)
            throw std::runtime_error("weighted MISF level exceeded the allocated hierarchy");
        misfSize[misfLevel] = mishSizeCurrLevel;
        misfLevel++;
    } while(itr);

    if(misfLevel > log_2_n)
        throw std::runtime_error("weighted MISF level exceeded the allocated hierarchy");
    misfLevel -= 1;
    while(misfLevel > 0 && misfSize[misfLevel] < numOfInitVert)
        misfLevel--;

    size_tt v = 0;
    if(misfSize[misfLevel] > numOfInitVert){
        while(v < numOfInitVert)
            vertDepth[mish[v++]] = misfLevel + 1;
        misfLevel++;
        misfSize[misfLevel] = numOfInitVert;
    }

    initMishHeight = misfLevel;
}

double DrawGraph::metric_me_init_v1(size_tt root)
{
    size_tt rootDepth = vertDepth[root];
    metricNbrs[root] = new std::vector<MetricNeighbor>[rootDepth + 1];
    metricNbrsDepth[root] = rootDepth + 1;

    size_tt bottomNbrsLayer = 0;
    double maxDist = 0.0;
    traverse_weighted_shortest_paths(
        root,
        std::numeric_limits<double>::infinity(),
        0,
        [&](size_tt overt, double distOvert){
            maxDist = std::max(maxDist, distOvert);
            size_tt limit = std::min(vertDepth[overt], rootDepth);
            for(size_tt i = bottomNbrsLayer; i <= limit; i++){
                if(metricNbrs[root][i].size() < nbr[i]){
                    metricNbrs[root][i].push_back(MetricNeighbor{overt, distOvert});
                } else {
                    bottomNbrsLayer = i + 1;
                }
            }
            return false;
        });

    return maxDist;
}

void DrawGraph::select_insertion_anchor_subset_weighted(std::vector<size_tt> &anchors,
                                                        std::vector<double> &anchorDist,
                                                        size_tt targetCount)
{
    if(targetCount == 0 || anchors.size() <= targetCount)
        return;

    Point<> candidateCentroid = pos[anchors[0]];
    candidateCentroid.set_to_zero();
    for(size_t i = 0; i < anchors.size(); i++)
        candidateCentroid += pos[anchors[i]];
    candidateCentroid /= (coord_t)anchors.size();

    if(insertionAnchorStrategy == INSERT_ANCHOR_STRATEGY_SPREAD_PREV &&
       dim == 2){
        const double pi = 3.14159265358979323846;
        struct PolarCandidate {
            size_t idx;
            double angle;
            double radius2;
            size_tt vert;
        };

        auto angular_gap = [pi](double a, double b){
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
        for(size_t i = 0; i < anchors.size(); i++){
            const Point<> &p = pos[anchors[i]];
            double dx = p.getX() - candidateCentroid.getX();
            double dy = p.getY() - candidateCentroid.getY();
            polar.push_back(PolarCandidate{
                i,
                std::atan2(dy, dx),
                dx * dx + dy * dy,
                anchors[i]
            });
        }

        std::vector<size_tt> selected;
        std::vector<double> selectedDist;
        std::vector<bool> used(anchors.size(), false);

        size_t firstChoice = 0;
        for(size_t i = 1; i < polar.size(); i++){
            if(polar[i].radius2 > polar[firstChoice].radius2 + 1e-12 ||
               (std::fabs(polar[i].radius2 - polar[firstChoice].radius2) <= 1e-12 &&
                polar[i].vert < polar[firstChoice].vert)){
                firstChoice = i;
            }
        }
        selected.push_back(anchors[polar[firstChoice].idx]);
        selectedDist.push_back(anchorDist[polar[firstChoice].idx]);
        used[polar[firstChoice].idx] = true;

        while(selected.size() < targetCount){
            size_t bestPolarIndex = polar.size();
            double bestAngleScore = -1.0;
            double bestSepScore = -1.0;
            for(size_t i = 0; i < polar.size(); i++){
                if(used[polar[i].idx])
                    continue;
                double minAngleGap = std::numeric_limits<double>::infinity();
                double minSep = std::numeric_limits<double>::infinity();
                for(size_t j = 0; j < selected.size(); j++){
                    size_tt chosenVert = selected[j];
                    const Point<> &chosen = pos[chosenVert];
                    double cdx = chosen.getX() - candidateCentroid.getX();
                    double cdy = chosen.getY() - candidateCentroid.getY();
                    double chosenAngle = std::atan2(cdy, cdx);
                    minAngleGap = std::min(minAngleGap,
                                           angular_gap(polar[i].angle, chosenAngle));
                    minSep = std::min(minSep,
                                      norm2(pos[anchors[polar[i].idx]], chosen));
                }
                if(minAngleGap > bestAngleScore + 1e-12 ||
                   (std::fabs(minAngleGap - bestAngleScore) <= 1e-12 &&
                    (minSep > bestSepScore + 1e-12 ||
                     (std::fabs(minSep - bestSepScore) <= 1e-12 &&
                      (bestPolarIndex >= polar.size() ||
                       polar[i].vert < polar[bestPolarIndex].vert))))){
                    bestAngleScore = minAngleGap;
                    bestSepScore = minSep;
                    bestPolarIndex = i;
                }
            }
            if(bestPolarIndex >= polar.size())
                break;
            used[polar[bestPolarIndex].idx] = true;
            selected.push_back(anchors[polar[bestPolarIndex].idx]);
            selectedDist.push_back(anchorDist[polar[bestPolarIndex].idx]);
        }

        if(selected.size() == targetCount){
            anchors.swap(selected);
            anchorDist.swap(selectedDist);
            return;
        }
    }

    if(insertionAnchorStrategy == INSERT_ANCHOR_STRATEGY_BALANCED_BAND){
        auto subsetObjective = [&](const std::vector<size_t> &subsetIndices){
            Point<> subsetCentroid = candidateCentroid;
            subsetCentroid.set_to_zero();
            for(size_t idx : subsetIndices)
                subsetCentroid += pos[anchors[idx]];
            subsetCentroid /= (coord_t)subsetIndices.size();

            double centroidError = norm2(subsetCentroid, candidateCentroid);
            double minSep = 0.0;
            if(subsetIndices.size() >= 2){
                minSep = std::numeric_limits<double>::infinity();
                for(size_t i = 0; i < subsetIndices.size(); i++){
                    for(size_t j = i + 1; j < subsetIndices.size(); j++){
                        double sep = norm2(pos[anchors[subsetIndices[i]]],
                                           pos[anchors[subsetIndices[j]]]);
                        if(sep < minSep)
                            minSep = sep;
                    }
                }
                if(!std::isfinite(minSep))
                    minSep = 0.0;
            }

            double scale = 0.0;
            for(size_t i = 0; i < anchors.size(); i++)
                scale = std::max(scale, norm2(pos[anchors[i]], candidateCentroid));
            if(scale <= 0.0)
                scale = 1.0;

            return centroidError / scale - 0.25 * (minSep / scale);
        };

        auto build_subset = [&](const std::vector<size_t> &bestIndices){
            std::vector<size_tt> selected;
            std::vector<double> selectedDist;
            selected.reserve(bestIndices.size());
            selectedDist.reserve(bestIndices.size());
            for(size_t idx : bestIndices){
                selected.push_back(anchors[idx]);
                selectedDist.push_back(anchorDist[idx]);
            }
            anchors.swap(selected);
            anchorDist.swap(selectedDist);
        };

        size_t n = anchors.size();
        size_t k = std::min<size_t>(targetCount, n);
        double combinationCount = 1.0;
        for(size_t i = 1; i <= k; i++){
            combinationCount *= static_cast<double>(n - k + i);
            combinationCount /= static_cast<double>(i);
        }

        if(combinationCount <= 50000.0){
            std::vector<size_t> current;
            std::vector<size_t> best;
            double bestObjective = std::numeric_limits<double>::infinity();

            std::function<void(size_t, size_t)> dfs =
                [&](size_t start, size_t remaining){
                    if(remaining == 0){
                        double objective = subsetObjective(current);
                        if(objective < bestObjective - 1e-12){
                            bestObjective = objective;
                            best = current;
                        } else if(std::fabs(objective - bestObjective) <= 1e-12){
                            if(best.empty() ||
                               anchors[current[0]] < anchors[best[0]]){
                                best = current;
                            }
                        }
                        return;
                    }
                    size_t limit = n - remaining;
                    for(size_t idx = start; idx <= limit; idx++){
                        current.push_back(idx);
                        dfs(idx + 1, remaining - 1);
                        current.pop_back();
                    }
                };
            dfs(0, k);
            if(!best.empty()){
                build_subset(best);
                return;
            }
        }

        std::vector<size_t> greedy;
        greedy.reserve(k);
        std::vector<bool> used(anchors.size(), false);
        while(greedy.size() < k){
            size_t bestIndex = anchors.size();
            double bestObjective = std::numeric_limits<double>::infinity();
            for(size_t i = 0; i < anchors.size(); i++){
                if(used[i])
                    continue;
                greedy.push_back(i);
                double objective = subsetObjective(greedy);
                greedy.pop_back();
                if(objective < bestObjective - 1e-12 ||
                   (std::fabs(objective - bestObjective) <= 1e-12 &&
                    (bestIndex >= anchors.size() || anchors[i] < anchors[bestIndex]))){
                    bestObjective = objective;
                    bestIndex = i;
                }
            }
            if(bestIndex >= anchors.size())
                break;
            used[bestIndex] = true;
            greedy.push_back(bestIndex);
        }
        if(!greedy.empty()){
            build_subset(greedy);
            return;
        }
    }

    std::vector<size_tt> selected;
    std::vector<double> selectedDist;
    std::vector<bool> used(anchors.size(), false);

    size_t firstChoice = 0;
    double firstScore = -1.0;
    for(size_t i = 0; i < anchors.size(); i++){
        double score = norm2(pos[anchors[i]], candidateCentroid);
        if(score > firstScore ||
           (std::fabs(score - firstScore) <= 1e-12 &&
            anchors[i] < anchors[firstChoice])){
            firstScore = score;
            firstChoice = i;
        }
    }
    selected.push_back(anchors[firstChoice]);
    selectedDist.push_back(anchorDist[firstChoice]);
    used[firstChoice] = true;

    while(selected.size() < targetCount){
        size_t bestIndex = anchors.size();
        double bestScore = -1.0;
        for(size_t i = 0; i < anchors.size(); i++){
            if(used[i])
                continue;
            double minSep = std::numeric_limits<double>::infinity();
            for(size_t j = 0; j < selected.size(); j++){
                double sep = norm2(pos[anchors[i]], pos[selected[j]]);
                if(sep < minSep)
                    minSep = sep;
            }
            if(minSep > bestScore ||
               (std::fabs(minSep - bestScore) <= 1e-12 &&
                (bestIndex >= anchors.size() || anchors[i] < anchors[bestIndex]))){
                bestScore = minSep;
                bestIndex = i;
            }
        }
        if(bestIndex >= anchors.size())
            break;
        used[bestIndex] = true;
        selected.push_back(anchors[bestIndex]);
        selectedDist.push_back(anchorDist[bestIndex]);
    }

    anchors.swap(selected);
    anchorDist.swap(selectedDist);
}

Point<> DrawGraph::initial_position_mode_weighted(const size_tt *closeVert,
                                                  const double *closeVertDist,
                                                  size_tt count,
                                                  size_tt placement_mode)
{
    if(placement_mode == LEVEL0_INSERT_LEAST_SQUARES)
        return initial_position_least_squares_weighted(closeVert, closeVertDist, count);
    if(placement_mode == PLACEMENT_CIRCLE && dim == 2)
        return initial_position_circle_weighted(closeVert, closeVertDist, count);
    return initial_position_barycenter(closeVert, count);
}

Point<> DrawGraph::initial_position_circle_weighted(const size_tt *closeVert,
                                                    const double *closeVertDist,
                                                    size_tt count)
{
    if(count < 3)
        return initial_position_barycenter(closeVert, count);

    auto dist2d = [](const Point<> &a, const Point<> &b){
        coord_t dx = a.getX() - b.getX();
        coord_t dy = a.getY() - b.getY();
        return std::sqrt(dx * dx + dy * dy);
    };

    auto add_intersections = [&](const Point<> &p1, coord_t r1,
                                 const Point<> &p2, coord_t r2,
                                 std::vector<Point<> > &out){
        coord_t dx = p2.getX() - p1.getX();
        coord_t dy = p2.getY() - p1.getY();
        coord_t d = std::sqrt(dx * dx + dy * dy);
        if(d <= 0)
            return;
        if(d > r1 + r2)
            return;
        if(d < std::abs(r1 - r2))
            return;

        coord_t a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
        coord_t h2 = r1 * r1 - a * a;
        if(h2 < 0)
            h2 = 0;
        coord_t h = std::sqrt(h2);
        coord_t xm = p1.getX() + a * dx / d;
        coord_t ym = p1.getY() + a * dy / d;
        coord_t rx = -dy * (h / d);
        coord_t ry = dx * (h / d);

        out.emplace_back(xm + rx, ym + ry, 0.0);
        if(h2 > 0)
            out.emplace_back(xm - rx, ym - ry, 0.0);
    };

    std::vector<Point<> > candidates;
    coord_t r0 = closeVertDist[0] * edge;
    coord_t r1 = closeVertDist[1] * edge;
    coord_t r2 = closeVertDist[2] * edge;

    add_intersections(pos[closeVert[0]], r0, pos[closeVert[1]], r1, candidates);
    add_intersections(pos[closeVert[1]], r1, pos[closeVert[2]], r2, candidates);
    add_intersections(pos[closeVert[0]], r0, pos[closeVert[2]], r2, candidates);

    if(candidates.size() < 3)
        return initial_position_barycenter(closeVert, count);

    coord_t best_score = std::numeric_limits<coord_t>::infinity();
    Point<> best;
    for(size_t i = 0; i < candidates.size(); i++){
        for(size_t j = i + 1; j < candidates.size(); j++){
            for(size_t k = j + 1; k < candidates.size(); k++){
                coord_t score = dist2d(candidates[i], candidates[j]) +
                                dist2d(candidates[i], candidates[k]) +
                                dist2d(candidates[j], candidates[k]);
                if(score < best_score){
                    best_score = score;
                    best = candidates[i] + candidates[j] + candidates[k];
                    best /= 3.0;
                }
            }
        }
    }

    if(best_score == std::numeric_limits<coord_t>::infinity())
        return initial_position_barycenter(closeVert, count);
    return best;
}

Point<> DrawGraph::initial_position_least_squares_weighted(const size_tt *closeVert,
                                                           const double *closeVertDist,
                                                           size_tt count)
{
    if(count <= dim)
        return initial_position_barycenter(closeVert, count);

    const Point<> &base = pos[closeVert[0]];
    coord_t base_coords[3] = {base.getX(), base.getY(), base.getZ()};
    coord_t d0 = closeVertDist[0] * edge;
    coord_t base_norm2 = 0.0;
    for(size_tt k = 0; k < dim; k++)
        base_norm2 += base_coords[k] * base_coords[k];

    double ata[3][3] = {{0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
    double atb[3] = {0.0, 0.0, 0.0};

    for(size_tt i = 1; i < count; i++){
        const Point<> &anchor = pos[closeVert[i]];
        coord_t anchor_coords[3] = {anchor.getX(), anchor.getY(), anchor.getZ()};
        coord_t di = closeVertDist[i] * edge;
        coord_t anchor_norm2 = 0.0;
        double row[3] = {0.0, 0.0, 0.0};
        for(size_tt k = 0; k < dim; k++){
            row[k] = 2.0 * (anchor_coords[k] - base_coords[k]);
            anchor_norm2 += anchor_coords[k] * anchor_coords[k];
        }
        double rhs = anchor_norm2 - base_norm2 - di * di + d0 * d0;
        for(size_tt rLocal = 0; rLocal < dim; rLocal++){
            atb[rLocal] += row[rLocal] * rhs;
            for(size_tt cLocal = 0; cLocal < dim; cLocal++)
                ata[rLocal][cLocal] += row[rLocal] * row[cLocal];
        }
    }

    double aug[3][4] = {{0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0}};
    for(size_tt rLocal = 0; rLocal < dim; rLocal++){
        for(size_tt cLocal = 0; cLocal < dim; cLocal++)
            aug[rLocal][cLocal] = ata[rLocal][cLocal];
        aug[rLocal][dim] = atb[rLocal];
    }

    const double tol = 1e-10;
    for(size_tt col = 0; col < dim; col++){
        size_tt pivot = col;
        for(size_tt row = col + 1; row < dim; row++){
            if(std::fabs(aug[row][col]) > std::fabs(aug[pivot][col]))
                pivot = row;
        }
        if(std::fabs(aug[pivot][col]) <= tol)
            return initial_position_barycenter(closeVert, count);
        if(pivot != col){
            for(size_tt k = col; k <= dim; k++)
                std::swap(aug[col][k], aug[pivot][k]);
        }
        double pivot_val = aug[col][col];
        for(size_tt k = col; k <= dim; k++)
            aug[col][k] /= pivot_val;
        for(size_tt row = 0; row < dim; row++){
            if(row == col)
                continue;
            double factor = aug[row][col];
            if(factor == 0.0)
                continue;
            for(size_tt k = col; k <= dim; k++)
                aug[row][k] -= factor * aug[col][k];
        }
    }

    coord_t sol[3] = {0.0, 0.0, 0.0};
    for(size_tt k = 0; k < dim; k++)
        sol[k] = aug[k][dim];
    return Point<>(sol[0], sol[1], dim > 2 ? sol[2] : 0.0);
}

void DrawGraph::metric_me_v1(size_tt root)
{
    bool level0Insertion = (misfLevel == 0);
    size_tt numOfCloseVert = level0Insertion ? level0AnchorCount : insertionAnchorCount;
    if(numOfCloseVert == 0)
        numOfCloseVert = 1;
    std::vector<size_tt> closeVert;
    std::vector<double> closeVertDist;
    closeVert.reserve(numOfCloseVert);
    closeVertDist.reserve(numOfCloseVert);
    double closeVertCutoffDist = 0.0;
    bool closeVertDone = false;
    size_tt insertionPlacementMode =
        level0Insertion ? level0InsertionMode : LEVEL0_INSERT_INHERIT;
    if(insertionPlacementMode == LEVEL0_INSERT_INHERIT)
        insertionPlacementMode = placementMode;
    size_tt localKkSteps = level0Insertion ? level0LocalKkSteps : 3;
    size_tt rootDepth = vertDepth[root];

    auto eligibleAnchor = [&](size_tt overt){
        if(vertDepth[overt] <= rootDepth)
            return false;
        if(insertionAnchorScope == INSERT_ANCHOR_SCOPE_PREV_MISF)
            return vertDepth[overt] == rootDepth + 1;
        return true;
    };

    auto finalizeInsertion = [&](size_tt anchorCount){
        if(anchorCount == 0){
            pos[root] = rand_Point();
            oldDisp[root].set_to_zero();
            oldDispNorm[root] = 0;
            record_insertion_trace(root,
                                   traceLevelIndex + 1,
                                   misfLevel,
                                   prevSize,
                                   misfSize[misfLevel],
                                   rootDepth,
                                   numOfCloseVert,
                                   0,
                                   insertionPlacementMode,
                                   localKkSteps,
                                   closeVert,
                                   closeVertDist,
                                   pos[root],
                                   pos[root],
                                   oldDisp[root],
                                   oldDisp[root],
                                   oldDispNorm[root],
                                   oldDispNorm[root]);
            return;
        }
        select_insertion_anchor_subset_weighted(closeVert, closeVertDist, numOfCloseVert);
        anchorCount = std::min<size_tt>(anchorCount, closeVert.size());
        pos[root] = initial_position_mode_weighted(closeVert.data(),
                                                   closeVertDist.data(),
                                                   anchorCount,
                                                   insertionPlacementMode);
        Point<> coordInitial = pos[root];

        oldDisp[root].set_to_zero();
        for(size_tt i = 0; i < anchorCount; i++)
            oldDisp[root] += oldDisp[closeVert[i]];
        oldDisp[root] /= (coord_t)anchorCount;
        oldDispNorm[root] = oldDisp[root].norm();
        Point<> oldDispInitial = oldDisp[root];
        coord_t oldDispNormInitial = oldDispNorm[root];

        for(size_tt itr = 0; itr < localKkSteps; itr++){
            KK_spring_weighted_local_v1(root,
                                        closeVert.data(),
                                        closeVertDist.data(),
                                        anchorCount);
            update_Local_Temp_v3(root, r, s);
            oldDisp[root] = disp[root];
            oldDispNorm[root] = dispNorm[root];
            disp[root] *= (coord_t)heat[root];
            if(dispNorm[root])
                disp[root] /= dispNorm[root];
            pos[root] += disp[root];
        }
        record_insertion_trace(root,
                               traceLevelIndex + 1,
                               misfLevel,
                               prevSize,
                               misfSize[misfLevel],
                               rootDepth,
                               numOfCloseVert,
                               anchorCount,
                               insertionPlacementMode,
                               localKkSteps,
                               closeVert,
                               closeVertDist,
                               coordInitial,
                               pos[root],
                               oldDispInitial,
                               oldDisp[root],
                               oldDispNormInitial,
                               oldDispNorm[root]);
    };

    metricNbrs[root] = new std::vector<MetricNeighbor>[rootDepth + 1];
    metricNbrsDepth[root] = rootDepth + 1;

    size_tt bottomNbrsLayer = 0;
    traverse_weighted_shortest_paths(
        root,
        std::numeric_limits<double>::infinity(),
        metricNeighborCap,
        [&](size_tt overt, double distOvert){
            size_tt limit = std::min(vertDepth[overt], rootDepth);
            for(size_tt i = bottomNbrsLayer; i <= limit; i++){
                if(metricNbrs[root][i].size() < nbr[i]){
                    metricNbrs[root][i].push_back(MetricNeighbor{overt, distOvert});
                } else {
                    bottomNbrsLayer = i + 1;
                }
            }

            if(!closeVertDone && eligibleAnchor(overt)){
                closeVert.push_back(overt);
                closeVertDist.push_back(distOvert);

                if(insertionAnchorStrategy == INSERT_ANCHOR_STRATEGY_FIRST){
                    if(closeVert.size() == numOfCloseVert){
                        closeVertDone = true;
                        finalizeInsertion(numOfCloseVert);
                    }
                } else if(closeVert.size() == numOfCloseVert){
                    closeVertCutoffDist = distOvert;
                }
            }

            if(!closeVertDone &&
               insertionAnchorStrategy != INSERT_ANCHOR_STRATEGY_FIRST &&
               closeVertCutoffDist > 0.0 &&
               distOvert > closeVertCutoffDist + kMetricTol){
                closeVertDone = true;
                finalizeInsertion(closeVert.size());
            }

            return closeVertDone && bottomNbrsLayer > rootDepth;
        });

    if(!closeVertDone && !closeVert.empty())
        finalizeInsertion(closeVert.size());
    else if(!closeVertDone)
        finalizeInsertion(0);
}

void DrawGraph::KK_spring_weighted_local_v1(const size_tt vert,
                                            const size_tt *closeVert,
                                            const double *closeVertDist,
                                            size_tt size)
{
    coord_t norm2Local;

    disp[vert].set_to_zero();
    for(size_tt i = 0; i < size; i++){
        if(closeVertDist[i] <= 0.0)
            continue;
        vect.set_to_zero();
        vect = pos[closeVert[i]] - pos[vert];
        norm2Local = vect.norm2();
        vect *= ((float)norm2Local /
                 (closeVertDist[i] * closeVertDist[i] * get_Edge2()) - 1);
        disp[vert] += vect;
    }
    dispNorm[vert] = disp[vert].norm();
}

void DrawGraph::KK_spring_weighted_v1(const size_tt vert,
                                      const std::vector<MetricNeighbor> &rootNbrsLayer,
                                      size_tt mishLayer)
{
    disp[vert].set_to_zero();
    std::unique_ptr<std::ostringstream> attractionEdges;
    if(collectRefinementStepDetails){
        attractionEdges.reset(new std::ostringstream());
        *attractionEdges << std::setprecision(17);
    }
    for(const MetricNeighbor &neighbor : rootNbrsLayer){
        if(neighbor.vert >= numOfVert || neighbor.vert == vert || neighbor.dist <= 0.0)
            continue;
        vect.set_to_zero();
        vect += pos[neighbor.vert];
        vect -= pos[vert];
        double norm2Local = (double)vect.fnorm2();
        std::vector<double> deltaFlat;
        coord_t desired = 0;
        coord_t desired2 = 0;
        double finalScale = 0.0;
        if(collectRefinementStepDetails){
            deltaFlat.resize(dim);
            deltaFlat[0] = vect.getX();
            if(dim > 1)
                deltaFlat[1] = vect.getY();
            if(dim > 2)
                deltaFlat[2] = vect.getZ();
            desired = edge * neighbor.dist;
            desired2 = desired * desired;
            finalScale =
                static_cast<double>(static_cast<float>(norm2Local / desired2 - 1.0));
        }
        vect *= (float)(norm2Local / (neighbor.dist * neighbor.dist * edge2) - 1.0);
        if(collectRefinementStepDetails){
            if(attractionEdges->tellp() > 0)
                *attractionEdges << ";";
            *attractionEdges << static_cast<int>(neighbor.vert) + 1 << ":"
                             << neighbor.dist;
            lastAttractionDisp[0] += vect.getX();
            if(dim > 1)
                lastAttractionDisp[1] += vect.getY();
            if(dim > 2)
                lastAttractionDisp[2] += vect.getZ();
            std::vector<double> stepFlat(dim);
            std::vector<double> cumulativeFlat(dim);
            stepFlat[0] = vect.getX();
            cumulativeFlat[0] = lastAttractionDisp[0];
            if(dim > 1){
                stepFlat[1] = vect.getY();
                cumulativeFlat[1] = lastAttractionDisp[1];
            }
            if(dim > 2){
                stepFlat[2] = vect.getZ();
                cumulativeFlat[2] = lastAttractionDisp[2];
            }
            lastAttractionTermNeighbors.push_back(static_cast<int>(neighbor.vert) + 1);
            lastAttractionTermWeights.push_back(neighbor.dist);
            lastAttractionTermNorm2.push_back(norm2Local);
            lastAttractionTermDesired.push_back(desired);
            lastAttractionTermDesired2.push_back(desired2);
            lastAttractionTermScale.push_back(finalScale);
            lastAttractionTermDelta.push_back(std::move(deltaFlat));
            lastAttractionTermStep.push_back(std::move(stepFlat));
            lastAttractionTermCumulative.push_back(std::move(cumulativeFlat));
        }
        disp[vert] += vect;
    }
    if(collectRefinementStepDetails)
        lastAttractionEdges = attractionEdges->str();

    if(mishLayer > 0 && activeVertCount > 1)
        add_coarse_global_repulsion(vert, activeVertCount);

    coord_t norm = disp[vert].fnorm();
    dispNorm[vert] = grip::detail::round_to_integer<unsigned long>(norm);

    if(dispNorm[vert]){
        disp[vert] *= edge / norm;
        dispNorm[vert] = disp[vert].norm();
    }
}

void DrawGraph::KK_spring_weighted_final_v1(const size_tt vert,
                                            const std::vector<MetricNeighbor> &rootNbrsLayer,
                                            size_tt mishLayer)
{
    (void)mishLayer;
    disp[vert].set_to_zero();
    std::unique_ptr<std::ostringstream> attractionEdges;
    if(collectRefinementStepDetails){
        attractionEdges.reset(new std::ostringstream());
        *attractionEdges << std::setprecision(17);
    }
    for(const MetricNeighbor &neighbor : rootNbrsLayer){
        if(neighbor.vert >= numOfVert || neighbor.vert == vert || neighbor.dist <= 0.0)
            continue;
        vect.set_to_zero();
        vect += pos[neighbor.vert];
        vect -= pos[vert];
        double norm2Local = (double)vect.fnorm2();
        std::vector<double> deltaFlat;
        coord_t desired = 0;
        coord_t desired2 = 0;
        double finalScale = 0.0;
        if(collectRefinementStepDetails){
            deltaFlat.resize(dim);
            deltaFlat[0] = vect.getX();
            if(dim > 1)
                deltaFlat[1] = vect.getY();
            if(dim > 2)
                deltaFlat[2] = vect.getZ();
            desired = edge * neighbor.dist;
            desired2 = desired * desired;
            finalScale =
                static_cast<double>(static_cast<float>(norm2Local / desired2 - 1.0));
        }
        vect *= (float)(norm2Local / (neighbor.dist * neighbor.dist * edge2) - 1.0);
        if(collectRefinementStepDetails){
            if(attractionEdges->tellp() > 0)
                *attractionEdges << ";";
            *attractionEdges << static_cast<int>(neighbor.vert) + 1 << ":"
                             << neighbor.dist;
            lastAttractionDisp[0] += vect.getX();
            if(dim > 1)
                lastAttractionDisp[1] += vect.getY();
            if(dim > 2)
                lastAttractionDisp[2] += vect.getZ();
            std::vector<double> stepFlat(dim);
            std::vector<double> cumulativeFlat(dim);
            stepFlat[0] = vect.getX();
            cumulativeFlat[0] = lastAttractionDisp[0];
            if(dim > 1){
                stepFlat[1] = vect.getY();
                cumulativeFlat[1] = lastAttractionDisp[1];
            }
            if(dim > 2){
                stepFlat[2] = vect.getZ();
                cumulativeFlat[2] = lastAttractionDisp[2];
            }
            lastAttractionTermNeighbors.push_back(static_cast<int>(neighbor.vert) + 1);
            lastAttractionTermWeights.push_back(neighbor.dist);
            lastAttractionTermNorm2.push_back(norm2Local);
            lastAttractionTermDesired.push_back(desired);
            lastAttractionTermDesired2.push_back(desired2);
            lastAttractionTermScale.push_back(finalScale);
            lastAttractionTermDelta.push_back(std::move(deltaFlat));
            lastAttractionTermStep.push_back(std::move(stepFlat));
            lastAttractionTermCumulative.push_back(std::move(cumulativeFlat));
        }
        disp[vert] += vect;
    }
    if(collectRefinementStepDetails)
        lastAttractionEdges = attractionEdges->str();

    if(activeVertCount > 1)
        add_active_global_repulsion(vert, activeVertCount, fedge2);

    coord_t norm = disp[vert].fnorm();
    dispNorm[vert] = grip::detail::round_to_integer<unsigned long>(norm);

    if(dispNorm[vert]){
        disp[vert] *= edge / norm;
        dispNorm[vert] = disp[vert].norm();
    }
}

void DrawGraph::FR_spring_weighted_v1(const size_tt vert,
                                      const std::vector<MetricNeighbor> &vertNbrs,
                                      size_tt misfLayer)
{
    size_tt overt;
    double norm2Local;

    disp[vert].set_to_zero();
    if(collectRefinementStepDetails){
        std::fill(lastAttractionDisp.begin(), lastAttractionDisp.end(), 0.0);
        std::fill(lastRepulsionDisp.begin(), lastRepulsionDisp.end(), 0.0);
        lastAttractionTermNeighbors.clear();
        lastAttractionTermWeights.clear();
        lastAttractionTermNorm2.clear();
        lastAttractionTermDesired.clear();
        lastAttractionTermDesired2.clear();
        lastAttractionTermScale.clear();
        lastAttractionTermDelta.clear();
        lastAttractionTermStep.clear();
        lastAttractionTermCumulative.clear();
    }
    std::unique_ptr<std::ostringstream> attractionEdges;
    std::unique_ptr<std::ostringstream> repulsionNeighbors;
    if(collectRefinementStepDetails){
        attractionEdges.reset(new std::ostringstream());
        repulsionNeighbors.reset(new std::ostringstream());
        *attractionEdges << std::setprecision(17);
        *repulsionNeighbors << std::setprecision(17);
    }

    size_tt degLocal = graph.adjList[0][vert];
    for(size_tt adjVert = 0; adjVert < degLocal; adjVert++){
        overt = graph.adjList[vert + 1][adjVert];
        vect.set_to_zero();
        vect += pos[overt];
        vect -= pos[vert];
        norm2Local = (double)vect.norm2();
        std::vector<double> deltaFlat;
        if(collectRefinementStepDetails){
            deltaFlat.resize(dim);
            deltaFlat[0] = vect.getX();
            if(dim > 1)
                deltaFlat[1] = vect.getY();
            if(dim > 2)
                deltaFlat[2] = vect.getZ();
        }
        vect *= (float)norm2Local;
        coord_t edgeWeight = graph.get_edge_weight(vert, adjVert);
        coord_t desired = edge * edgeWeight;
        coord_t desired2 = desired * desired;
        vect /= desired2;
        double finalScale = 0.0;
        if(collectRefinementStepDetails){
            finalScale =
                static_cast<double>(static_cast<float>(norm2Local)) / desired2;
            if(attractionEdges->tellp() > 0)
                *attractionEdges << ";";
            *attractionEdges << static_cast<int>(overt) + 1 << ":"
                             << edgeWeight;
            lastAttractionDisp[0] += vect.getX();
            if(dim > 1)
                lastAttractionDisp[1] += vect.getY();
            if(dim > 2)
                lastAttractionDisp[2] += vect.getZ();
            std::vector<double> stepFlat(dim);
            std::vector<double> cumulativeFlat(dim);
            stepFlat[0] = vect.getX();
            cumulativeFlat[0] = lastAttractionDisp[0];
            if(dim > 1){
                stepFlat[1] = vect.getY();
                cumulativeFlat[1] = lastAttractionDisp[1];
            }
            if(dim > 2){
                stepFlat[2] = vect.getZ();
                cumulativeFlat[2] = lastAttractionDisp[2];
            }
            lastAttractionTermNeighbors.push_back(static_cast<int>(overt) + 1);
            lastAttractionTermWeights.push_back(edgeWeight);
            lastAttractionTermNorm2.push_back(norm2Local);
            lastAttractionTermDesired.push_back(desired);
            lastAttractionTermDesired2.push_back(desired2);
            lastAttractionTermScale.push_back(finalScale);
            lastAttractionTermDelta.push_back(std::move(deltaFlat));
            lastAttractionTermStep.push_back(std::move(stepFlat));
            lastAttractionTermCumulative.push_back(std::move(cumulativeFlat));
        }
        disp[vert] += vect;
    }

    for(const MetricNeighbor &neighbor : vertNbrs){
        if(neighbor.vert >= numOfVert || neighbor.vert == vert || neighbor.dist <= 0.0)
            continue;
        vect.set_to_zero();
        vect += pos[vert];
        vect -= pos[neighbor.vert];
        norm2Local = (double)vect.fnorm2();
        if(!norm2Local)
            continue;
        vect *= (float)(fedge2 / norm2Local);
        if(collectRefinementStepDetails){
            if(repulsionNeighbors->tellp() > 0)
                *repulsionNeighbors << ";";
            *repulsionNeighbors << static_cast<int>(neighbor.vert) + 1 << ":"
                                << neighbor.dist;
            lastRepulsionDisp[0] += vect.getX();
            if(dim > 1)
                lastRepulsionDisp[1] += vect.getY();
            if(dim > 2)
                lastRepulsionDisp[2] += vect.getZ();
        }
        disp[vert] += vect;
    }

    if(!misfLayer)
        add_final_anchor_force(vert);

    coord_t norm = disp[vert].fnorm();
    dispNorm[vert] = grip::detail::round_to_integer<unsigned long>(norm);

    if(dispNorm[vert]){
        disp[vert] *= edge / norm;
        dispNorm[vert] = disp[vert].norm();
    }
    if(collectRefinementStepDetails){
        lastAttractionEdges = attractionEdges->str();
        lastRepulsionNeighbors = repulsionNeighbors->str();
    }
}

void DrawGraph::mish_engine_weighted()
{
    bool firstRound = true;
    size_tt csize = numOfInitVert;
    size_tt ctr = 0;
    size_tt trace_round_in_level = 0;
    bool loop = true;

    while(loop && prevSize != csize){
        if(displayPar)
            loop = false;

        if(firstRound){
            firstRound = false;
            rounds = initRounds;
            prevMishSize = 0;
            currentRoundInLevel = 0;

            for(size_tt i = 0; i < csize; i++){
                pos[mish[i]] = rand_Point();
                center += pos[mish[i]];
            }
            baricenter = center / (coord_t)numOfInitVert;
            for(size_tt i = 0; i < csize; i++)
                pos[mish[i]] -= baricenter;
            trace_round_in_level = 0;
            if(!misfLevel && csize == numOfVert){
                for(size_tt i = 0; i < csize; i++)
                    finalAnchorPos[mish[i]] = pos[mish[i]];
                finalAnchorReady = true;
            }
            trace_begin_level(csize);
        } else if(ctr == rounds){
            if(!createList)
                lgkk_refine_level(csize, misfLevel, trace_round_in_level);
            ctr = 0;
            prevSize = csize;
            prevMishSize = misfSize[misfLevel];

            if(prevMishSize < numOfInitVert)
                prevMishSize = numOfInitVert;
            if(misfSize[misfLevel] != numOfVert)
                csize = misfSize[--misfLevel];
            else if(prevSize != numOfVert)
                csize = numOfVert;
            else
                break;

            rounds = sched3(csize,
                            0, initRounds,
                            numOfVert, finalRounds);

            for(size_tt i = 0; i < prevSize; i++)
                heat[mish[i]] = tinit;

            for(size_tt i = prevSize; i < csize; i++)
                metric_me_v1(mish[i]);
            currentRoundInLevel = 0;
            if(!misfLevel && csize == numOfVert){
                for(size_tt i = 0; i < csize; i++)
                    finalAnchorPos[mish[i]] = pos[mish[i]];
                finalAnchorReady = true;
            }
            trace_round_in_level = 0;
            trace_begin_level(csize);
        }

        if(ctr++ < rounds){
            activeVertCount = csize;
            currentRoundInLevel = ctr;
            bool recordSteps =
                should_record_refinement_step(traceLevelIndex, misfLevel, currentRoundInLevel);
            std::vector<size_t> stepRows;
            if(recordSteps)
                stepRows.assign(csize, static_cast<size_t>(-1));
            collectRefinementStepDetails = recordSteps;
            for(size_tt i = 0; i < csize; i++){
                Point<> coordBefore;
                coord_t heatBefore = 0;
                float oldCosBefore = 0;
                coord_t oldDispNormBefore = 0;
                if(recordSteps){
                    coordBefore = pos[mish[i]];
                    heatBefore = heat[mish[i]];
                    oldCosBefore = old_cos[mish[i]];
                    oldDispNormBefore = oldDispNorm[mish[i]];
                    std::fill(lastAttractionDisp.begin(), lastAttractionDisp.end(), 0.0);
                    std::fill(lastRepulsionDisp.begin(), lastRepulsionDisp.end(), 0.0);
                    lastAttractionEdges.clear();
                    lastRepulsionNeighbors.clear();
                    lastAttractionTermNeighbors.clear();
                    lastAttractionTermWeights.clear();
                    lastAttractionTermNorm2.clear();
                    lastAttractionTermDesired.clear();
                    lastAttractionTermDesired2.clear();
                    lastAttractionTermScale.clear();
                    lastAttractionTermDelta.clear();
                    lastAttractionTermStep.clear();
                    lastAttractionTermCumulative.clear();
                }
                if(!misfLevel && finalStageMode == FINAL_STAGE_KK_REPULSE)
                    KK_spring_weighted_final_v1(mish[i], metricNbrs[mish[i]][misfLevel], misfLevel);
                else if(!misfLevel)
                    FR_spring_weighted_v1(mish[i], metricNbrs[mish[i]][misfLevel], misfLevel);
                else
                    KK_spring_weighted_v1(mish[i], metricNbrs[mish[i]][misfLevel], misfLevel);

                Point<> preTempDisp;
                coord_t preTempDispNorm = 0;
                if(recordSteps){
                    preTempDisp = disp[mish[i]];
                    preTempDispNorm = dispNorm[mish[i]];
                }
                update_Local_Temp_v2(mish[i]);
                coord_t heatAfter = 0;
                float oldCosAfter = 0;
                if(recordSteps){
                    heatAfter = heat[mish[i]];
                    oldCosAfter = old_cos[mish[i]];
                }
                oldDisp[mish[i]] = disp[mish[i]];
                oldDispNorm[mish[i]] = dispNorm[mish[i]];
                disp[mish[i]] *= (coord_t)heat[mish[i]];

                if(dispNorm[mish[i]])
                    disp[mish[i]] /= dispNorm[mish[i]];
                if(!misfLevel &&
                   finalStageMode == FINAL_STAGE_FR &&
                   currentRoundInLevel > 1 &&
                   finalMoveScaleAfterFirst < 1)
                    disp[mish[i]] *= finalMoveScaleAfterFirst;
                if(recordSteps){
                    stepRows[i] = record_refinement_step_pre(
                        mish[i],
                        i + 1,
                        csize,
                        traceLevelIndex,
                        misfLevel,
                        currentRoundInLevel,
                        coordBefore,
                        preTempDisp,
                        lastAttractionDisp,
                        lastRepulsionDisp,
                        disp[mish[i]],
                        heatBefore,
                        heatAfter,
                        oldCosBefore,
                        oldCosAfter,
                        oldDispNormBefore,
                        preTempDispNorm,
                        lastAttractionEdges,
                        lastRepulsionNeighbors
                    );
                }
            }
            for(size_tt i = 0; i < csize; i++){
                pos[mish[i]] += disp[mish[i]];
                if(recordSteps && stepRows[i] != static_cast<size_t>(-1))
                    record_refinement_step_after(stepRows[i], mish[i]);
            }
            collectRefinementStepDetails = false;
            trace_round_in_level = ctr;
            trace_after_round(csize, trace_round_in_level);
        }
    }

    trace_finalize(csize, trace_round_in_level);

    if(loop && listSwitch){
        createList = true;
        listSwitch = false;
    } else {
        createList = false;
    }
}
