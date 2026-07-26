// Shared support routines for the GRIP layout engine

#include <algorithm>
#include <limits>
#include <queue>
#include <unordered_set>

#include "DrawGraph.h"


//****************************************************************
//
//       bfs_me_init_v2
//
//       Initialize retained graph-distance neighborhoods with BFS.
//
//****************************************************************
size_tt DrawGraph::bfs_me_init_v2(size_tt root)
{
    size_tt bottomNbrsLayer = 0; // the index i of the first free nbrs[root][i]
    size_tt vert;
    size_tt currentDepth = 0;
    FastQueue<size_tt> vertDepthQueue(4*numOfVert);// an array based queue
    marked[root] = root;
    vertDepthQueue.enqueue(root);
    vertDepthQueue.enqueue(currentDepth);
    
    // memory allocation for rootNbrs and nbrCounter
    size_tt *nbrCounter = new size_tt[vertDepth[root] + 1];
    nbrs[root] = new size_tt*[vertDepth[root] + 1];
    nbrsDepth[root] = vertDepth[root] + 1;
    for(size_tt i=0; i <= vertDepth[root]; i++){
        // 1 slot for vertex and 1 (the next one)
        // for its distance from the root
        nbrs[root][i] = new size_tt[2*nbr[i]]();
        nbrCounter[i] = 0;
    }

    do{
        vert = vertDepthQueue.dequeue();
        currentDepth = vertDepthQueue.dequeue() + 1;
        
        size_tt deg = graph.adjList[0][vert];
        for (size_tt adjVert = 0; adjVert < deg; adjVert++){
            size_tt overt = graph.adjList[vert+1][adjVert];

            if ( marked[overt] != (int)root ){
                marked[overt] = root;
                vertDepthQueue.enqueue(overt);
                vertDepthQueue.enqueue(currentDepth);

                for(int i = bottomNbrsLayer;
                    i <= std::min(vertDepth[overt], vertDepth[root]); i++){
                    if(nbrCounter[i] < 2*nbr[i]){
                        nbrs[root][i][nbrCounter[i]++] = overt;
                        nbrs[root][i][nbrCounter[i]++] = currentDepth;
                    } else
                        bottomNbrsLayer = i+1;
                }                
            }            
        }
    } while ( bottomNbrsLayer <= vertDepth[root] &&
              !vertDepthQueue.is_empty() );

    delete [] nbrCounter;
    return currentDepth;
}

//****************************************************************
//
//       bfs_me_v4
//
//       Place a vertex and retain its graph-distance neighborhood with BFS.
//
//****************************************************************
void DrawGraph::bfs_me_v4(size_tt root)
{
    size_tt bottomNbrsLayer = 0; // the index i of the first free nbrs[root][i]
    size_tt vert;
    size_tt currentDepth = 0;
    FastQueue<size_tt> vertDepthQueue(4*numOfVert);// an array based queue

    bool level0Insertion = (misfLevel == 0);
    size_tt numOfCloseVert = level0Insertion ? level0AnchorCount : insertionAnchorCount;
    if(numOfCloseVert == 0)
        numOfCloseVert = 1;
    std::vector<size_tt> closeVert;
    std::vector<size_tt> closeVertDist;
    closeVert.reserve(numOfCloseVert);
    closeVertDist.reserve(numOfCloseVert);
    size_tt closeVertCutoffDepth = 0;
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
        if(anchorCount == 0)
            return;
        select_insertion_anchor_subset(closeVert, closeVertDist, numOfCloseVert);
        anchorCount = std::min<size_tt>(anchorCount, closeVert.size());
        pos[root] = initial_position_mode(closeVert.data(),
                                          closeVertDist.data(),
                                          anchorCount,
                                          insertionPlacementMode);

        oldDisp[root].set_to_zero();
        for(size_tt i = 0; i < anchorCount; i++)
            oldDisp[root] += oldDisp[closeVert[i]];
        oldDisp[root] /= (coord_t)anchorCount;
        oldDispNorm[root] = oldDisp[root].norm();

        for(size_tt itr = 0; itr < localKkSteps; itr++){
            KK_spring_local(root, closeVert.data(),
                            closeVertDist.data(), anchorCount);
            update_Local_Temp_v3(root, r, s);
            oldDisp[root] = disp[root];
            oldDispNorm[root] = dispNorm[root];
            disp[root] *= (coord_t)heat[root];
            if(dispNorm[root])
                disp[root] /= dispNorm[root];
            pos[root] += disp[root];
        }
    };

    marked[root] = root;
    vertDepthQueue.enqueue(root);
    vertDepthQueue.enqueue(currentDepth);
    
    // memory allocation for rootNbrs
    nbrs[root] = new size_tt*[vertDepth[root]+1];
    nbrsDepth[root] = vertDepth[root] + 1;
    size_tt *nbrCounter = new size_tt[vertDepth[root]+1];
    for(size_tt i=0; i <= vertDepth[root]; i++){
        // 1 slot for vertex and 1 (the next one)
        // for its distance from the root
        nbrs[root][i] = new size_tt[2*nbr[i]]();
        nbrCounter[i] = 0;
    }
    
    do{
        vert = vertDepthQueue.dequeue();
        currentDepth = vertDepthQueue.dequeue() + 1;
        
        for (size_tt adjVert = 0; adjVert < graph.adjList[0][vert]; adjVert++){
            size_tt overt = graph.adjList[vert+1][adjVert];
                
            if ( marked[overt] != (int)root ){
                marked[overt] = root;
                vertDepthQueue.enqueue(overt);
                vertDepthQueue.enqueue(currentDepth);
                
                for(int i = bottomNbrsLayer;
                    i <= std::min(vertDepth[overt], vertDepth[root]); i++)
                    if(nbrCounter[i] < 2*nbr[i]){
                        nbrs[root][i][nbrCounter[i]++] = overt;
                        nbrs[root][i][nbrCounter[i]++] = currentDepth;
                    } else
                        bottomNbrsLayer = i+1;

                if(!closeVertDone && eligibleAnchor(overt)){
                    closeVertDist.push_back(currentDepth);
                    closeVert.push_back(overt);

                    if(insertionAnchorStrategy == INSERT_ANCHOR_STRATEGY_FIRST){
                        if(closeVert.size() == numOfCloseVert){
                            closeVertDone = true;
                            finalizeInsertion(numOfCloseVert);
                        }
                    } else if(closeVert.size() == numOfCloseVert){
                        closeVertCutoffDepth = currentDepth;
                    }
                }//end of if( !closeVertDone ...
            }            
        }
        if(!closeVertDone &&
           insertionAnchorStrategy != INSERT_ANCHOR_STRATEGY_FIRST &&
           closeVertCutoffDepth > 0 &&
           currentDepth > closeVertCutoffDepth){
            closeVertDone = true;
            finalizeInsertion(closeVert.size());
        }
    } while ( (!closeVertDone || bottomNbrsLayer <= vertDepth[root]) &&
              !vertDepthQueue.is_empty( ) );

    if(!closeVertDone && !closeVert.empty())
        finalizeInsertion(closeVert.size());
    delete [] nbrCounter;
}

void DrawGraph::compute_active_shortest_paths(size_tt sourceIndex,
                                              size_tt activeCount,
                                              std::vector<double> &dist,
                                              std::vector<int> *parent)
{
    const double inf = std::numeric_limits<double>::infinity();
    const double tol = 1e-10;
    dist.assign(activeCount, inf);
    if(parent)
        parent->assign(activeCount, -1);
    if(sourceIndex >= activeCount)
        return;

    dist[sourceIndex] = 0.0;

    if(!graph.has_weights()){
        std::queue<size_tt> q;
        std::vector<char> inQueue(activeCount, 0);
        q.push(sourceIndex);
        inQueue[sourceIndex] = 1;

        while(!q.empty()){
            size_tt currentIndex = q.front();
            q.pop();
            inQueue[currentIndex] = 0;
            size_tt currentVert = mish[currentIndex];
            double currentDist = dist[currentIndex];
            size_tt deg = graph.adjList[0][currentVert];
            for(size_tt adjVert = 0; adjVert < deg; adjVert++){
                size_tt overt = graph.adjList[currentVert + 1][adjVert];
                int overtIndex = lgkkActiveIndex[overt];
                if(overtIndex < 0 || static_cast<size_tt>(overtIndex) >= activeCount)
                    continue;
                double alt = currentDist + 1.0;
                double best = dist[overtIndex];
                bool improve = alt + tol < best;
                bool equal = std::isfinite(best) && std::fabs(alt - best) <= tol;
                if(improve){
                    dist[overtIndex] = alt;
                    if(parent)
                        (*parent)[overtIndex] = static_cast<int>(currentIndex);
                    if(!inQueue[overtIndex]){
                        q.push(static_cast<size_tt>(overtIndex));
                        inQueue[overtIndex] = 1;
                    }
                } else if(equal && parent){
                    int currentParent = (*parent)[overtIndex];
                    if(currentParent < 0 ||
                       mish[currentIndex] < mish[static_cast<size_tt>(currentParent)]){
                        (*parent)[overtIndex] = static_cast<int>(currentIndex);
                    }
                }
            }
        }
        return;
    }

    struct QueueNode {
        double dist;
        size_tt vert;
        size_tt index;
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
    pq.push(QueueNode{0.0, mish[sourceIndex], sourceIndex});

    while(!pq.empty()){
        double currentDist = pq.top().dist;
        size_tt currentIndex = pq.top().index;
        pq.pop();
        if(currentDist > dist[currentIndex] + tol)
            continue;
        size_tt currentVert = mish[currentIndex];
        size_tt deg = graph.adjList[0][currentVert];
        for(size_tt adjVert = 0; adjVert < deg; adjVert++){
            size_tt overt = graph.adjList[currentVert + 1][adjVert];
            int overtIndex = lgkkActiveIndex[overt];
            if(overtIndex < 0 || static_cast<size_tt>(overtIndex) >= activeCount)
                continue;
            double alt = currentDist + graph.get_edge_weight(currentVert, adjVert);
            double best = dist[overtIndex];
            double scale = std::max(1.0,
                                    std::max(std::fabs(alt),
                                             std::isfinite(best) ? std::fabs(best) : 0.0));
            bool improve = !std::isfinite(best) || alt + tol * scale < best;
            bool equal = std::isfinite(best) &&
                std::fabs(alt - best) <= tol * scale;
            if(improve){
                dist[overtIndex] = alt;
                if(parent)
                    (*parent)[overtIndex] = static_cast<int>(currentIndex);
                pq.push(QueueNode{alt, overt, static_cast<size_tt>(overtIndex)});
            } else if(equal && parent){
                int currentParent = (*parent)[overtIndex];
                if(currentParent < 0 ||
                   mish[currentIndex] < mish[static_cast<size_tt>(currentParent)]){
                    (*parent)[overtIndex] = static_cast<int>(currentIndex);
                }
            }
        }
    }
}

std::vector<size_tt> DrawGraph::lgkk_choose_local_neighbors(size_tt sourceIndex,
                                                            size_tt activeCount) const
{
    struct Candidate {
        double dist;
        size_tt index;
        size_tt vert;
    };
    std::vector<Candidate> candidates;
    if(lgkkLocalNbrs == 0 || sourceIndex >= activeCount)
        return std::vector<size_tt>();

    const double *row = lgkkDistanceMatrix.data() + sourceIndex * activeCount;
    candidates.reserve(activeCount > 0 ? activeCount - 1 : 0);
    for(size_tt idx = 0; idx < activeCount; idx++){
        if(idx == sourceIndex || !std::isfinite(row[idx]))
            continue;
        candidates.push_back(Candidate{row[idx], idx, mish[idx]});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate &lhs, const Candidate &rhs){
                  if(lhs.dist != rhs.dist)
                      return lhs.dist < rhs.dist;
                  return lhs.vert < rhs.vert;
              });

    std::vector<size_tt> out;
    out.reserve(std::min<size_tt>(lgkkLocalNbrs, candidates.size()));
    for(size_t i = 0; i < candidates.size() && out.size() < lgkkLocalNbrs; i++){
        out.push_back(candidates[i].index);
    }
    return out;
}

std::vector<size_tt> DrawGraph::lgkk_choose_landmarks(size_tt sourceIndex,
                                                      size_tt activeCount) const
{
    if(lgkkLandmarkCount == 0 || sourceIndex >= activeCount)
        return std::vector<size_tt>();

    std::vector<size_tt> candidates;
    candidates.reserve(activeCount > 0 ? activeCount - 1 : 0);
    const double *sourceRow = lgkkDistanceMatrix.data() + sourceIndex * activeCount;
    for(size_tt idx = 0; idx < activeCount; idx++){
        if(idx == sourceIndex || !std::isfinite(sourceRow[idx]))
            continue;
        candidates.push_back(idx);
    }
    if(candidates.empty())
        return std::vector<size_tt>();

    std::vector<size_tt> selected;
    std::vector<double> coverage;
    coverage.reserve(candidates.size());
    for(size_t i = 0; i < candidates.size(); i++)
        coverage.push_back(sourceRow[candidates[i]]);

    size_t maxCount = std::min<size_t>(lgkkLandmarkCount, candidates.size());
    for(size_t step = 0; step < maxCount; step++){
        size_t choicePos = candidates.size();
        double choiceScore = -1.0;
        for(size_t i = 0; i < candidates.size(); i++){
            double score = (step == 0) ? sourceRow[candidates[i]] : coverage[i];
            if(score > choiceScore ||
               (std::fabs(score - choiceScore) <= 1e-12 &&
                (choicePos >= candidates.size() ||
                 mish[candidates[i]] < mish[candidates[choicePos]]))){
                choiceScore = score;
                choicePos = i;
            }
        }
        if(choicePos >= candidates.size())
            break;
        size_tt choice = candidates[choicePos];
        selected.push_back(choice);
        candidates.erase(candidates.begin() + static_cast<long>(choicePos));
        coverage.erase(coverage.begin() + static_cast<long>(choicePos));
        if(candidates.empty())
            break;
        const double *choiceRow = lgkkDistanceMatrix.data() + choice * activeCount;
        for(size_t i = 0; i < candidates.size(); i++)
            coverage[i] = std::min(coverage[i], choiceRow[candidates[i]]);
    }
    return selected;
}

void DrawGraph::build_lgkk_level_cache(size_tt activeCount,
                                       size_tt mishLayer)
{
    clear_lgkk_level_cache();
    if(!should_run_multiscale_lgkk(activeCount, mishLayer))
        return;

    lgkkActiveIndex.assign(numOfVert, -1);
    for(size_tt i = 0; i < activeCount; i++)
        lgkkActiveIndex[mish[i]] = static_cast<int>(i);

    lgkkDistanceMatrix.assign(activeCount * activeCount,
                              std::numeric_limits<double>::infinity());
    std::vector<double> dist;
    for(size_tt sourceIndex = 0; sourceIndex < activeCount; sourceIndex++){
        compute_active_shortest_paths(sourceIndex, activeCount, dist, nullptr);
        std::copy(dist.begin(),
                  dist.end(),
                  lgkkDistanceMatrix.begin() + sourceIndex * activeCount);
    }

    std::vector<std::vector<size_tt>> selectedTargets(activeCount);
    for(size_tt sourceIndex = 0; sourceIndex < activeCount; sourceIndex++){
        std::vector<size_tt> local = lgkk_choose_local_neighbors(sourceIndex, activeCount);
        std::vector<size_tt> landmarks = lgkk_choose_landmarks(sourceIndex, activeCount);
        local.insert(local.end(), landmarks.begin(), landmarks.end());
        std::sort(local.begin(), local.end());
        local.erase(std::unique(local.begin(), local.end()), local.end());
        selectedTargets[sourceIndex] = std::move(local);
    }

    std::unordered_set<uint64_t> seenPairs;
    std::vector<int> parent;
    lgkkPairs.reserve(activeCount * std::max<size_tt>(1, lgkkLocalNbrs + lgkkLandmarkCount));
    for(size_tt sourceIndex = 0; sourceIndex < activeCount; sourceIndex++){
        if(selectedTargets[sourceIndex].empty())
            continue;
        compute_active_shortest_paths(sourceIndex, activeCount, dist, &parent);
        for(size_t targetPos = 0; targetPos < selectedTargets[sourceIndex].size(); targetPos++){
            size_tt targetIndex = selectedTargets[sourceIndex][targetPos];
            size_tt sourceVert = mish[sourceIndex];
            size_tt targetVert = mish[targetIndex];
            size_tt minVert = std::min(sourceVert, targetVert);
            size_tt maxVert = std::max(sourceVert, targetVert);
            uint64_t key = (static_cast<uint64_t>(minVert) << 32) |
                           static_cast<uint64_t>(maxVert);
            if(seenPairs.find(key) != seenPairs.end())
                continue;
            if(!std::isfinite(dist[targetIndex]))
                continue;

            std::vector<size_tt> pathVertices;
            int currentIndex = static_cast<int>(targetIndex);
            while(currentIndex >= 0 && static_cast<size_tt>(currentIndex) != sourceIndex){
                pathVertices.push_back(mish[static_cast<size_tt>(currentIndex)]);
                currentIndex = parent[static_cast<size_tt>(currentIndex)];
            }
            if(currentIndex < 0)
                continue;
            pathVertices.push_back(sourceVert);
            std::reverse(pathVertices.begin(), pathVertices.end());
            if(pathVertices.size() < 2)
                continue;

            LgkkPairCache pair;
            pair.source = sourceVert;
            pair.target = targetVert;
            pair.graphDistance = static_cast<coord_t>(dist[targetIndex]);
            pair.pathEdges.reserve(pathVertices.size() - 1);
            for(size_t edgeIndex = 1; edgeIndex < pathVertices.size(); edgeIndex++){
                pair.pathEdges.push_back(
                    LgkkPathEdge{pathVertices[edgeIndex - 1], pathVertices[edgeIndex]}
                );
            }
            lgkkPairs.push_back(std::move(pair));
            seenPairs.insert(key);
        }
    }

    double numerator = 0.0;
    double denominator = 0.0;
    const double eps2 = 1e-16;
    for(size_t pairIndex = 0; pairIndex < lgkkPairs.size(); pairIndex++){
        const LgkkPairCache &pair = lgkkPairs[pairIndex];
        double h = 0.0;
        for(size_t edgeIndex = 0; edgeIndex < pair.pathEdges.size(); edgeIndex++){
            const LgkkPathEdge &edgeRef = pair.pathEdges[edgeIndex];
            vect.set_to_zero();
            vect += pos[edgeRef.u];
            vect -= pos[edgeRef.v];
            h += std::sqrt(vect.fnorm2() + eps2);
        }
        double g = std::max<double>(pair.graphDistance, 1e-8);
        double kk = 1.0 / (g * g);
        numerator += kk * g * h;
        denominator += kk * g * g;
    }

    lgkkCacheScaleL0 = denominator > 0.0
        ? static_cast<coord_t>(numerator / denominator)
        : 1.0;
    lgkkCacheActiveCount = activeCount;
    lgkkCacheMisfLevel = static_cast<int>(mishLayer);
}

void DrawGraph::lgkk_refine_level(size_tt activeCount,
                                  size_tt mishLayer,
                                  size_tt &traceRoundInLevel)
{
    if(!should_run_multiscale_lgkk(activeCount, mishLayer))
        return;
    size_tt roundBudget = lgkk_round_budget_for_layer(mishLayer);
    if(roundBudget == 0)
        return;
    build_lgkk_level_cache(activeCount, mishLayer);
    if(lgkkPairs.empty())
        return;

    const double eps2 = 1e-16;
    const double initialStep = 1.0;
    const double stepShrink = 0.5;
    const double armijo = 1e-4;
    const double gradTol2 = 1e-16;
    const double minStep = 1e-8;
    const double distanceFloor = 1e-8;

    struct LgkkState {
        double energy;
        double gradNorm2;
        std::vector<Point<>> gradient;
    };

    std::vector<Point<>> activePos(activeCount);
    for(size_tt i = 0; i < activeCount; i++)
        activePos[i] = pos[mish[i]];

    auto evaluate_state = [&](const std::vector<Point<>> &coords){
        LgkkState state;
        state.energy = 0.0;
        state.gradNorm2 = 0.0;
        state.gradient.assign(activeCount, Point<>());
        for(size_tt i = 0; i < activeCount; i++)
            state.gradient[i].set_to_zero();

        for(size_t pairIndex = 0; pairIndex < lgkkPairs.size(); pairIndex++){
            const LgkkPairCache &pair = lgkkPairs[pairIndex];
            double g = std::max<double>(pair.graphDistance, distanceFloor);
            double kk = 1.0 / (g * g);
            double target = lgkkCacheScaleL0 * pair.graphDistance;

            std::vector<Point<>> edgeDiffs;
            std::vector<double> edgeLens;
            edgeDiffs.reserve(pair.pathEdges.size());
            edgeLens.reserve(pair.pathEdges.size());

            double h = 0.0;
            for(size_t edgeIndex = 0; edgeIndex < pair.pathEdges.size(); edgeIndex++){
                const LgkkPathEdge &edgeRef = pair.pathEdges[edgeIndex];
                int uIndex = lgkkActiveIndex[edgeRef.u];
                int vIndex = lgkkActiveIndex[edgeRef.v];
                if(uIndex < 0 || vIndex < 0)
                    continue;
                Point<> diff = coords[static_cast<size_t>(uIndex)] -
                               coords[static_cast<size_t>(vIndex)];
                double len = std::sqrt(diff.fnorm2() + eps2);
                edgeDiffs.push_back(diff);
                edgeLens.push_back(len);
                h += len;
            }
            if(edgeDiffs.empty())
                continue;

            double resid = h - target;
            double coeff = kk * resid;
            state.energy += 0.5 * kk * resid * resid;
            for(size_t edgeIndex = 0; edgeIndex < pair.pathEdges.size(); edgeIndex++){
                const LgkkPathEdge &edgeRef = pair.pathEdges[edgeIndex];
                int uIndex = lgkkActiveIndex[edgeRef.u];
                int vIndex = lgkkActiveIndex[edgeRef.v];
                if(uIndex < 0 || vIndex < 0)
                    continue;
                if(edgeLens[edgeIndex] <= 0.0)
                    continue;
                Point<> stepVec = edgeDiffs[edgeIndex] * (coeff / edgeLens[edgeIndex]);
                state.gradient[static_cast<size_t>(uIndex)] += stepVec;
                state.gradient[static_cast<size_t>(vIndex)] -= stepVec;
            }
        }

        for(size_tt i = 0; i < activeCount; i++)
            state.gradNorm2 += state.gradient[i].fnorm2();

        return state;
    };

    std::vector<Point<>> acceptedMove(activeCount);
    for(size_tt i = 0; i < activeCount; i++)
        acceptedMove[i].set_to_zero();

    LgkkState state = evaluate_state(activePos);
    for(size_tt roundIndex = 1; roundIndex <= roundBudget; roundIndex++){
        if(!std::isfinite(state.energy) || state.gradNorm2 <= gradTol2)
            break;

        double step = initialStep;
        bool accepted = false;
        std::vector<Point<>> proposal(activeCount);
        LgkkState candidate = state;

        while(std::isfinite(step) && step >= minStep){
            for(size_tt i = 0; i < activeCount; i++)
                proposal[i] = activePos[i] - state.gradient[i] * step;

            candidate = evaluate_state(proposal);
            double targetEnergy = state.energy - armijo * step * state.gradNorm2;
            if(std::isfinite(candidate.energy) && candidate.energy <= targetEnergy){
                accepted = true;
                break;
            }
            step *= stepShrink;
        }

        if(!accepted)
            break;

        for(size_tt i = 0; i < activeCount; i++){
            acceptedMove[i] = proposal[i] - activePos[i];
            activePos[i] = proposal[i];
            size_tt vert = mish[i];
            pos[vert] = activePos[i];
            disp[vert] = acceptedMove[i];
            oldDisp[vert] = acceptedMove[i];
            dispNorm[vert] = ROUND_L(acceptedMove[i].fnorm());
            oldDispNorm[vert] = dispNorm[vert];
        }

        state = candidate;
        currentRoundInLevel = rounds + roundIndex;
        traceRoundInLevel = currentRoundInLevel;
        if(traceMode == TRACE_ROUND &&
           (traceEvery <= 1 || roundIndex % traceEvery == 0)){
            capture_trace_snapshot("lgkk", activeCount, traceRoundInLevel);
        }
    }

    if(traceMode == TRACE_LEVEL)
        capture_trace_snapshot("lgkk", activeCount, traceRoundInLevel);
}


//**************************************************************
//
//    add_active_global_repulsion()
//
//    optional active-set repulsion over currently active vertices
//
//**************************************************************
void DrawGraph::add_active_global_repulsion(const size_tt vert,
                                            size_tt activeCount,
                                            coord_t repulsionScale)
{
    if(repulsionScale <= 0 || activeCount <= 1 || coarseRepulsionSample == 0)
        return;

    size_tt population = activeCount - 1;
    if(activeCount <= coarseRepulsionExactBelow ||
       coarseRepulsionSample >= population){
        add_active_global_repulsion_exact(vert, activeCount, repulsionScale);
        return;
    }

    add_active_global_repulsion_sampled(
        vert,
        activeCount,
        std::min(coarseRepulsionSample, population),
        repulsionScale
    );
}

void DrawGraph::add_active_global_repulsion_exact(const size_tt vert,
                                                  size_tt activeCount,
                                                  coord_t repulsionScale)
{
    for(size_tt i = 0; i < activeCount; i++){
        size_tt overt = mish[i];
        if(overt == vert)
            continue;
        vect.set_to_zero();
        vect += pos[vert];
        vect -= pos[overt];
        double norm2 = (double)vect.fnorm2();
        if(!norm2)
            continue;
        vect *= (float)(repulsionScale / norm2);
        if(collectRefinementStepDetails){
            if(lastRepulsionNeighbors.size() > 0)
                lastRepulsionNeighbors += ";";
            lastRepulsionNeighbors += std::to_string(static_cast<int>(overt) + 1);
            lastRepulsionDisp[0] += vect.getX();
            if(dim > 1)
                lastRepulsionDisp[1] += vect.getY();
            if(dim > 2)
                lastRepulsionDisp[2] += vect.getZ();
        }
        disp[vert] += vect;
    }
}

void DrawGraph::add_active_global_repulsion_sampled(const size_tt vert,
                                                    size_tt activeCount,
                                                    size_tt sampleCount,
                                                    coord_t repulsionScale)
{
    if(sampleCount == 0)
        return;

    std::vector<size_tt> sampled;
    sampled.reserve(sampleCount);
    while(sampled.size() < sampleCount){
        size_tt overt = mish[graph.fast_Rand() % activeCount];
        if(overt == vert)
            continue;
        if(std::find(sampled.begin(), sampled.end(), overt) != sampled.end())
            continue;
        sampled.push_back(overt);
    }

    double scale = (double)(activeCount - 1) / (double)sampleCount;
    for(size_tt overt : sampled){
        vect.set_to_zero();
        vect += pos[vert];
        vect -= pos[overt];
        double norm2 = (double)vect.fnorm2();
        if(!norm2)
            continue;
        vect *= (float)((repulsionScale * scale) / norm2);
        if(collectRefinementStepDetails){
            if(lastRepulsionNeighbors.size() > 0)
                lastRepulsionNeighbors += ";";
            lastRepulsionNeighbors += std::to_string(static_cast<int>(overt) + 1);
            lastRepulsionDisp[0] += vect.getX();
            if(dim > 1)
                lastRepulsionDisp[1] += vect.getY();
            if(dim > 2)
                lastRepulsionDisp[2] += vect.getZ();
        }
        disp[vert] += vect;
    }
}

void DrawGraph::add_coarse_global_repulsion(const size_tt vert,
                                            size_tt activeCount)
{
    add_active_global_repulsion(vert, activeCount, coarseFedge2);
}

void DrawGraph::add_coarse_global_repulsion_exact(const size_tt vert,
                                                  size_tt activeCount)
{
    add_active_global_repulsion_exact(vert, activeCount, coarseFedge2);
}

void DrawGraph::add_coarse_global_repulsion_sampled(const size_tt vert,
                                                    size_tt activeCount,
                                                    size_tt sampleCount)
{
    add_active_global_repulsion_sampled(vert,
                                        activeCount,
                                        sampleCount,
                                        coarseFedge2);
}

void DrawGraph::add_final_anchor_force(const size_tt vert)
{
    if(finalAnchorFactor <= 0 || !finalAnchorReady)
        return;

    vect.set_to_zero();
    vect += finalAnchorPos[vert];
    vect -= pos[vert];
    vect *= (float)(finalAnchorFactor / edge);
    disp[vert] += vect;
}

//**************************************************************
//
//    KK_spring_v4()
//
//    when misfLevel = 0 we use adjacent vertices for the force
//    calculation, otherwise use the whole bfsQNQueue
//    destroying it when preserveBFS is false.
//
//**************************************************************
void DrawGraph::KK_spring_v4(const size_tt vert,
                             size_tt *vertNbrs,
                             size_tt misfLayer)
{
    size_tt overt; // other vertex
    double dist2;// square of the graph theoretic dist between vert and overt
    double norm2; // Square of the Euclidean distance between vert and overt.
    
    size_tt *ptr;
    ptr = vertNbrs;

    disp[vert].set_to_zero();
    for(size_tt i = 0; i < 2*nbr[misfLayer]; i += 2){
        overt = *ptr++;
        dist2 = (double)(*ptr) * (*ptr);
        ptr++;
        if(!dist2 || overt >= numOfVert || overt == vert)
            continue;
        vect.set_to_zero();
        vect += pos[overt];
        vect -= pos[vert];
        norm2 = (double)vect.fnorm2();
        vect *= (float)(norm2/(dist2 * edge2) - 1);
        disp[vert] += vect;
    }

    if(misfLayer > 0 && activeVertCount > 1)
        add_coarse_global_repulsion(vert, activeVertCount);
    
    coord_t norm = disp[vert].fnorm();
    dispNorm[vert] = ROUND_L(norm);
    
    if(dispNorm[vert]){
        disp[vert] *= edge/norm;    
        dispNorm[vert] = disp[vert].norm();
    }
}

//**************************************************************
//
//    KK_spring_final()
//
//    finest-level local KK refinement with explicit active-set repulsion
//
//**************************************************************
void DrawGraph::KK_spring_final(const size_tt vert,
                                size_tt *vertNbrs,
                                size_tt misfLayer)
{
    size_tt overt;
    double dist2;
    double norm2;

    size_tt *ptr = vertNbrs;

    disp[vert].set_to_zero();
    for(size_tt i = 0; i < 2*nbr[misfLayer]; i += 2){
        overt = *ptr++;
        dist2 = (double)(*ptr) * (*ptr);
        ptr++;
        if(!dist2 || overt >= numOfVert || overt == vert)
            continue;
        vect.set_to_zero();
        vect += pos[overt];
        vect -= pos[vert];
        norm2 = (double)vect.fnorm2();
        vect *= (float)(norm2/(dist2 * edge2) - 1);
        disp[vert] += vect;
    }

    if(activeVertCount > 1)
        add_active_global_repulsion(vert, activeVertCount, fedge2);

    coord_t norm = disp[vert].fnorm();
    dispNorm[vert] = ROUND_L(norm);

    if(dispNorm[vert]){
        disp[vert] *= edge/norm;
        dispNorm[vert] = disp[vert].norm();
    }
}

//**************************************************************
//
//    KK_spring_local()
//
//    Compute a local spring displacement for a newly inserted vertex
//    using its selected anchors and their graph distances.
//
//**************************************************************
void DrawGraph::KK_spring_local(const size_tt vert,
                                size_tt *closeVert,
                                size_tt *closeVertDist,
                                size_tt size)
{
    coord_t norm2;   // its norm squared
    disp[vert].set_to_zero();
    for(size_tt i = 0; i < size; i++){
        vect.set_to_zero();
        vect = pos[closeVert[i]] - pos[vert];
        norm2 = vect.norm2();
        vect *=
            ((float)norm2/(closeVertDist[i]*closeVertDist[i]*get_Edge2())-1);
        disp[vert] += vect;
    }
    dispNorm[vert] = disp[vert].norm();
}

//**************************************************************
//
//	update_Local_Temp_v2()
//
//      updating heat[vert]
//      looking at the ratio oldNorm/newNorm and multiplying
//      old heat by this value and cos of the angle
//
//**************************************************************
void DrawGraph::update_Local_Temp_v2( size_tt vert )
{
    coord_t temp = heat[vert];
    coord_t normOldDisp = oldDispNorm[vert];
    coord_t normNewDisp = dispNorm[vert];
    
    if( normOldDisp != 0 && normNewDisp != 0 ){
        coord_t scalProd = disp[vert] * oldDisp[vert];
        coord_t cos = scalProd/(normOldDisp * normNewDisp);
        if( old_cos[vert] * cos > 0 )
            temp += (coord_t)(temp * s * cos * r);
        else
            temp += (coord_t)(temp * cos * r);
        
        old_cos[vert] = cos;
        heat[vert] = temp;
    }
}

//**************************************************************
//
//	update_Local_Temp_v3()
//
//      updating heat[vert]
//      looking at the ratio oldNorm/newNorm and multiplying
//      old heat by this value and cos of the angle
//
//**************************************************************
void DrawGraph::update_Local_Temp_v3( size_tt vert, coord_t r, coord_t s)
{
    coord_t temp = heat[vert];
    coord_t normOldDisp = oldDispNorm[vert];
    coord_t normNewDisp = dispNorm[vert];
    
    if( normOldDisp != 0 && normNewDisp != 0 ){
        coord_t scalProd = disp[vert] * oldDisp[vert];
        coord_t cos = scalProd/(normOldDisp * normNewDisp);
        if( old_cos[vert] * cos > 0 )
            temp += (coord_t)(temp * s * cos * r);
        else
            temp += (coord_t)(temp * cos * r);
        
        old_cos[vert] = cos;
        heat[vert] = temp;
    }
}

//**************************************************************
//
//    FR_spring_v2()
//
//    force calculation
//
//**************************************************************
void DrawGraph::FR_spring_v2(const size_tt vert,
                             size_tt *vertNbrs,
                             size_tt misfLayer)
{
    size_tt overt; // other vertex
    double norm2; // Square of the Euclidean distance between vert and overt.
    
    size_tt *ptr;
    ptr = vertNbrs;

    disp[vert].set_to_zero();

        // attractive force calculation
        size_tt deg = graph.adjList[0][vert];
        for (size_tt adjVert = 0; adjVert < deg; adjVert++){
            overt = graph.adjList[vert+1][adjVert];
            vect.set_to_zero();
            vect += pos[overt];
            vect -= pos[vert];
            norm2 = (double)vect.fnorm2();
            vect *= (float)vect.norm2();
            coord_t desired2 = edge2;
            if(graph.has_weights()){
                coord_t w = graph.get_edge_weight(vert, adjVert);
                coord_t desired = edge * w;
                desired2 = desired * desired;
            }
            vect /= desired2;
            disp[vert] += vect;
        }

    // repulsive force calculation
    size_tt locNbr = 2*nbr[misfLayer];
    for(size_tt i = 0; i < locNbr; i += 2){
        overt = *ptr++;
        size_tt graphDist = *ptr++;
        if(!graphDist || overt >= numOfVert || overt == vert)
            continue;
        vect.set_to_zero();
        vect += pos[vert];
        vect -= pos[overt];
        norm2 = (double)vect.fnorm2();
        if(!norm2)
            continue;
        vect *= (float)(fedge2/norm2);
        disp[vert] += vect;
    }
    coord_t norm = disp[vert].fnorm();
    dispNorm[vert] = ROUND_L(norm);

    if(dispNorm[vert]){
        disp[vert] *= edge/norm;    
        dispNorm[vert] = disp[vert].norm();
    }
}
