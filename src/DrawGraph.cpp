// FILE: DrawGraph.cpp - member function definitions for DrawGraph class.
//   
// April 1, 2000 - moved all changes info to the file changes
//
#include "DrawGraph.h"
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>
#include <limits>

//**************************************************************
//
//      class constructor
//
//**************************************************************
DrawGraph::DrawGraph(const Graph &_graph,
                     size_tt _dim       ,
                     size_tt _rounds    ,
                     size_tt _finalRounds,
                     size_tt _tinit_factor,
                     size_tt _numOfInitVert,
                     size_tt _numOfNbrs ,
                     coord_t _r           ,  //parameters of update_Local_Temp_v3()
                     coord_t _s           ,
                     coord_t _repulsionFactor,
                     size_tt _placementMode,
                     bool _displayPar   ,
                     size_tt _finalStageMode,
                     coord_t _coarseRepulsionFactor,
                     size_tt _coarseRepulsionSample,
                     size_tt _coarseRepulsionExactBelow,
                     coord_t _finalAnchorFactor,
                     coord_t _finalMoveScaleAfterFirst,
                     size_tt _insertionAnchorCount,
                     size_tt _insertionAnchorScope,
                     size_tt _insertionAnchorStrategy,
                     size_tt _level0InsertionMode,
                     size_tt _level0AnchorCount,
                     size_tt _level0LocalKkSteps,
                     size_tt _lgkkMultiscaleRounds,
                     size_tt _lgkkLocalNbrs,
                     size_tt _lgkkLandmarkCount,
                     size_tt _lgkkScope,
                     size_tt _lgkkActiveLimit)
: createList(false),
  graph(_graph),
  dim(_dim),
  rounds(_rounds),
  initRounds(_rounds),
  finalRounds(_finalRounds),
  numOfVert(_graph.get_numOfVert()),
  edge(32),
  edge2(edge * edge),
  tinit(edge / _tinit_factor),
  SMALL_DIST(32),  // Initialize with the desired value
  SMALL_DIST2(2 * SMALL_DIST + 1),  // Initialize based on SMALL_DIST
  numOfInitVert(std::min(_numOfInitVert, _graph.get_numOfVert())),
  numOfNbrs(_numOfNbrs),
  r(_r), s(_s),
  listSwitch(true),
  displayPar(_displayPar),
  placementMode((_placementMode == PLACEMENT_CIRCLE)
                    ? PLACEMENT_CIRCLE
                    : PLACEMENT_BARYCENTER),
  finalStageMode((_finalStageMode == FINAL_STAGE_KK_REPULSE)
                     ? FINAL_STAGE_KK_REPULSE
                     : FINAL_STAGE_FR),
  fedge2(_repulsionFactor * 0.05 * edge2),
  coarseFedge2(_coarseRepulsionFactor * 0.05 * edge2),
  coarseRepulsionSample(_coarseRepulsionSample),
  coarseRepulsionExactBelow(_coarseRepulsionExactBelow),
  finalAnchorFactor(_finalAnchorFactor),
  finalMoveScaleAfterFirst(_finalMoveScaleAfterFirst),
  insertionAnchorCount(std::max<size_tt>(1, _insertionAnchorCount)),
  insertionAnchorScope((_insertionAnchorScope == INSERT_ANCHOR_SCOPE_PREV_MISF)
                           ? INSERT_ANCHOR_SCOPE_PREV_MISF
                           : INSERT_ANCHOR_SCOPE_ANY_HIGHER),
  insertionAnchorStrategy((_insertionAnchorStrategy == INSERT_ANCHOR_STRATEGY_DISTANCE_BAND)
                              ? INSERT_ANCHOR_STRATEGY_DISTANCE_BAND
                              : (_insertionAnchorStrategy == INSERT_ANCHOR_STRATEGY_BALANCED_BAND)
                                    ? INSERT_ANCHOR_STRATEGY_BALANCED_BAND
                                    : INSERT_ANCHOR_STRATEGY_FIRST),
  level0InsertionMode((_level0InsertionMode == LEVEL0_INSERT_BARYCENTER)
                          ? LEVEL0_INSERT_BARYCENTER
                          : (_level0InsertionMode == LEVEL0_INSERT_LEAST_SQUARES)
                                ? LEVEL0_INSERT_LEAST_SQUARES
                                : LEVEL0_INSERT_INHERIT),
  level0AnchorCount(std::max<size_tt>(1, _level0AnchorCount)),
  level0LocalKkSteps(_level0LocalKkSteps),
  lgkkMultiscaleRounds(_lgkkMultiscaleRounds),
  lgkkLocalNbrs(_lgkkLocalNbrs),
  lgkkLandmarkCount(_lgkkLandmarkCount),
  lgkkScope((_lgkkScope == LGKK_SCOPE_COARSE)
                ? LGKK_SCOPE_COARSE
                : LGKK_SCOPE_ALL),
  lgkkActiveLimit(std::max<size_tt>(1, _lgkkActiveLimit)),
  activeVertCount(0),
  currentRoundInLevel(0),
  finalAnchorReady(false),
  traceMode(TRACE_NONE),
  traceEvery(1),
  traceLevelIndex(0),
  finalAnchorPos(_graph.get_numOfVert()),
  lgkkCacheActiveCount(0),
  lgkkCacheMisfLevel(-1),
  lgkkCacheScaleL0(1.0)
{
#define DEBUG 0
    if( dim != 2 && dim != 3 )
        throw std::runtime_error("only 2D and 3D layouts are supported");
#if DEBUG 
    debug("Entering constructor\n");
#endif

    pos          = new Point<>[numOfVert];
    disp         = new Point<>[numOfVert];
    oldDisp      = new Point<>[numOfVert];

    dispNorm     = new coord_t[numOfVert];
    oldDispNorm  = new coord_t[numOfVert];
    old_cos      = new float[numOfVert];
    heat         = new coord_t[numOfVert];
    deg          = new size_tt[numOfVert];
        
    for(size_tt vert=0; vert < numOfVert; vert++){
        disp[ vert ].set_to_zero();
        oldDisp[ vert ].set_to_zero();
        dispNorm[ vert ] = 0;
        oldDispNorm[ vert ] = 0;
        heat[ vert ] = tinit;
        old_cos[vert] = 1;
        deg[vert] = graph.adjList[0][vert];
    }

    numOfInitVert = std::min(numOfInitVert, numOfVert);

    //
    // GRIP engine state
    //

    AvgDeg = 0;  // USED ONLY FOR STATISTICS
    maxCxty = 0;
    initCxty = 10000;
    smallLevel = 0;
    prevSize = 0;
        
    //computing Avg(deg(G)), maxCxty
    for(size_tt vert = 0; vert < numOfVert; vert++)
        AvgDeg += deg[vert];
        
    maxCxty = (unsigned long)AvgDeg;
    if(maxCxty < initCxty)
        maxCxty = initCxty;
    AvgDeg /= (float)numOfVert;
        
    //debug("maxCxty="<< maxCxty<< ", AvgDeg="<< AvgDeg);

    // mish is the "maximal independent set hierarchy"
    mish = new size_tt[numOfVert]; // max ind set hierarchy
    inv  = new size_tt[numOfVert]; // "inverse" of mish
    for(size_tt i=0; i<numOfVert; i++)
        inv[i] = i, mish[i] = i;

    //
    //    CREATING MISF
    //
    // initializing marked array, which will be used in BFSs
    marked = new int[numOfVert];
    for(int i = 0; i < numOfVert; i++) 
        marked[i] = -1;
    
    log_2_n = ilog(numOfVert)+2;
//        debug("log_2_n="<< log_2_n-2);
    misfSize  = new size_tt[log_2_n];
    vertDepth = new size_tt[numOfVert];
    for(size_tt i = 0; i < log_2_n; i++)
        misfSize[i] = 0;
    for(size_tt i = 0; i < numOfVert; i++)
        vertDepth[i] = 0;

    if(numOfVert == numOfInitVert){
        misfSize[0] = numOfVert;
        misfLevel = 0;
    } else 
        create_misf();
    
    // resetting marked array values for BFSs in the main part of
    // the program - MishEngine(s).cpp
    for(int i = 0; i < numOfVert; i++){
//        debug("marked["<<i<<"]="<<marked[i]);
        marked[i] = -1;
    }
    
    // Computing smallLevel i.e. a level so that for each
    // misf level l >= smallLevel
    // misfSize[l] * misfSize[l] <= initCxty
    size_tt itr = 0;
    while( itr < log_2_n && misfSize[itr] ){
        if( (double)misfSize[itr] * misfSize[itr] - initCxty <= 0){
            smallLevel = itr;
            break;
        }
        itr++;
    }
        
    //debug("smallLevel="<<smallLevel);

    // numOfNbrs directly controls the retained local graph-distance
    // neighborhood size used by the refinement routines.
    nbr       = new size_tt[log_2_n];
    itr = 0;
    while( itr < log_2_n && misfSize[itr] ){
        size_tt available = misfSize[itr] > 0 ? misfSize[itr] - 1 : 0;
        nbr[itr] = std::min(numOfNbrs, available);
        itr++;
    }

#if 0
    std::cout << "<DrawGraph.cpp>["<<__LINE__<<"] nbr:       ";
    itr = 0;
    while( itr < log_2_n && misfSize[itr] ){
        std::cout << nbr[itr++] << ' ';
    }
    std::cout << "\n";
        
    std::cout << "<DrawGraph.cpp>["<<__LINE__<<"] misfSize: ";
    itr = 0;
    while( itr < log_2_n && misfSize[itr] ){
        std::cout << misfSize[itr++] << ' ';
    }
    std::cout << "\n";
#endif

//          std::cout << "<DrawGraph.cpp>[169] mish: ";
//          for(size_tt v = 0; v < min(40,(int)numOfVert); v++)
//              std::cout << mish[v] << ' ';
//          std::cout << endl;
        
//          std::cout << "<DrawGraph.cpp>[174] vertDepth: ";
//          for(size_tt v = 0; v < numOfVert; v++)
//              std::cout << vertDepth[v]<< ' ';
//          std::cout << endl;
        
            
    nbrs = new size_tt**[numOfVert];
    nbrsDepth = new size_tt[numOfVert];
    for(size_tt i = 0; i < numOfVert; i++){
        nbrs[i] = nullptr;
        nbrsDepth[i] = 0;
    }
    diam = 0;
    size_tt height = 0;
    for(size_tt i = 0; i < numOfInitVert; i++){
        height = bfs_me_init_v2(mish[i]);
        if( diam < height )
            diam = height;
    }
    //debug("diam="<<diam);
        
    boxSize = (coord_t)(edge * .7 * diam);
    box2Size = 2 * boxSize + 1;
        
#if DEBUG    
    debug("Leaving constructor\n");
#endif
}

void DrawGraph::configure_trace(size_tt mode, size_tt every)
{
    traceMode = mode;
    traceEvery = std::max<size_tt>(1, every);
    traceLevelIndex = 0;
    traceFrames.clear();
    tracePhases.clear();
    traceLevelIndices.clear();
    traceMisfLevels.clear();
    traceRounds.clear();
    traceActiveCounts.clear();
}

void DrawGraph::capture_trace_snapshot(const char *phase,
                                       size_tt activeCount,
                                       size_tt roundInLevel)
{
    if(traceMode == TRACE_NONE)
        return;

    activeCount = std::min(activeCount, numOfVert);

    std::vector<double> frame(numOfVert * dim,
                              std::numeric_limits<double>::quiet_NaN());
    for(size_tt i = 0; i < activeCount; i++){
        size_tt vert = mish[i];
        frame[vert] = pos[vert].getX();
        if(dim > 1)
            frame[vert + numOfVert] = pos[vert].getY();
        if(dim > 2)
            frame[vert + 2 * numOfVert] = pos[vert].getZ();
    }

    traceFrames.push_back(std::move(frame));
    tracePhases.emplace_back(phase);
    traceLevelIndices.push_back(static_cast<int>(traceLevelIndex));
    traceMisfLevels.push_back(static_cast<int>(misfLevel));
    traceRounds.push_back(static_cast<int>(roundInLevel));
    traceActiveCounts.push_back(static_cast<int>(activeCount));
}

void DrawGraph::trace_begin_level(size_tt activeCount)
{
    if(traceMode == TRACE_NONE)
        return;

    traceLevelIndex++;
    if(traceMode == TRACE_ROUND ||
       traceLevelIndex == 1 ||
       ((traceLevelIndex - 1) % traceEvery == 0)){
        capture_trace_snapshot(traceLevelIndex == 1 ? "init" : "level_start",
                               activeCount,
                               0);
    }
}

void DrawGraph::trace_after_round(size_tt activeCount, size_tt roundInLevel)
{
    if(traceMode != TRACE_ROUND)
        return;
    if(traceEvery > 1 && roundInLevel % traceEvery != 0)
        return;
    capture_trace_snapshot("round", activeCount, roundInLevel);
}

void DrawGraph::trace_finalize(size_tt activeCount, size_tt roundInLevel)
{
    if(traceMode == TRACE_NONE)
        return;
    capture_trace_snapshot("final", activeCount, roundInLevel);
}

//**************************************************************
//
//     destructor
//
//**************************************************************
DrawGraph::~DrawGraph(){
    delete [] dispNorm;
    delete [] oldDispNorm;
    delete [] old_cos;
    delete [] heat;

    delete [] pos;
    delete [] disp;
    delete [] oldDisp;
    delete [] deg;
    delete [] mish;
    delete [] inv;
    delete [] misfSize;

    //
    // MORE PROPERLY the first double for loop should be
    // executed in KK_spring when unused parts of nbrs[i]
    // are deleted
    //
    if(nbrs){
        for(size_tt i = 0; i < numOfVert; i++){
            if(nbrs[i]){
                for(size_tt j = 0; j < nbrsDepth[i]; j++)
                    delete [] nbrs[i][j];
                delete [] nbrs[i];
            }
        }
        delete [] nbrs;
    }
    delete [] nbrsDepth;

    delete [] vertDepth;
    delete [] nbr;
    delete [] marked;
}


//**************************************************************
//
//	dist()
//
//**************************************************************
coord_t DrawGraph::dist(const Point<> & p, const Point<> & q)
{
  return (coord_t)sqrt(norm2(p,q));
}

//**************************************************************
//
//	rand_Point()
//
//	random point in the box centered at
//      the origin of the coordinate system with the edge
//      2 * SIZE.
//
//**************************************************************
Point<> DrawGraph::rand_Point()
{
    if( dim == 2 )
        return
            Point<>((coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize,
                    (coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize, 0);
    else
        return
            Point<>((coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize,
                    (coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize,
//                   0);
                    (coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize);
}

Point<> DrawGraph::initial_position(const size_tt *closeVert,
                                    const size_tt *closeVertDist,
                                    size_tt count)
{
    return initial_position_mode(closeVert, closeVertDist, count, placementMode);
}

Point<> DrawGraph::initial_position_mode(const size_tt *closeVert,
                                         const size_tt *closeVertDist,
                                         size_tt count,
                                         size_tt placement_mode)
{
    if( placement_mode == LEVEL0_INSERT_LEAST_SQUARES )
        return initial_position_least_squares(closeVert, closeVertDist, count);
    if( placement_mode == PLACEMENT_CIRCLE && dim == 2 )
        return initial_position_circle(closeVert, closeVertDist, count);
    return initial_position_barycenter(closeVert, count);
}

Point<> DrawGraph::initial_position_barycenter(const size_tt *closeVert,
                                               size_tt count)
{
    if(count == 0)
        return rand_Point();
    Point<> result = pos[closeVert[0]];
    for(size_tt i = 1; i < count; i++)
        result += pos[closeVert[i]];
    result /= (coord_t)count;
    return result;
}

Point<> DrawGraph::initial_position_circle(const size_tt *closeVert,
                                           const size_tt *closeVertDist,
                                           size_tt count)
{
    if(count < 3)
        return initial_position_barycenter(closeVert, count);

    auto dist2d = [](const Point<> &a, const Point<> &b){
        coord_t dx = a.getX() - b.getX();
        coord_t dy = a.getY() - b.getY();
        return std::sqrt(dx*dx + dy*dy);
    };

    auto add_intersections = [&](const Point<> &p1, coord_t r1,
                                 const Point<> &p2, coord_t r2,
                                 std::vector<Point<> > &out){
        coord_t dx = p2.getX() - p1.getX();
        coord_t dy = p2.getY() - p1.getY();
        coord_t d = std::sqrt(dx*dx + dy*dy);
        if(d <= 0)
            return;
        if(d > r1 + r2)
            return;
        if(d < std::abs(r1 - r2))
            return;

        coord_t a = (r1*r1 - r2*r2 + d*d) / (2.0 * d);
        coord_t h2 = r1*r1 - a*a;
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

Point<> DrawGraph::initial_position_least_squares(const size_tt *closeVert,
                                                  const size_tt *closeVertDist,
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
        for(size_tt r = 0; r < dim; r++){
            atb[r] += row[r] * rhs;
            for(size_tt c = 0; c < dim; c++)
                ata[r][c] += row[r] * row[c];
        }
    }

    double aug[3][4] = {{0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0}};
    for(size_tt r = 0; r < dim; r++){
        for(size_tt c = 0; c < dim; c++)
            aug[r][c] = ata[r][c];
        aug[r][dim] = atb[r];
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

void DrawGraph::select_insertion_anchor_subset(std::vector<size_tt> &anchors,
                                               std::vector<size_tt> &anchorDist,
                                               size_tt targetCount)
{
    if(targetCount == 0 || anchors.size() <= targetCount)
        return;

    Point<> candidateCentroid = pos[anchors[0]];
    candidateCentroid.set_to_zero();
    for(size_t i = 0; i < anchors.size(); i++)
        candidateCentroid += pos[anchors[i]];
    candidateCentroid /= (coord_t)anchors.size();

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
            std::vector<size_tt> selectedDist;
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
    std::vector<size_tt> selectedDist;
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

bool DrawGraph::should_run_multiscale_lgkk(size_tt activeCount,
                                           size_tt mishLayer) const
{
    if(lgkkMultiscaleRounds == 0)
        return false;
    if(lgkkLocalNbrs == 0 && lgkkLandmarkCount == 0)
        return false;
    if(activeCount < 2 || activeCount > lgkkActiveLimit)
        return false;
    if(lgkkScope == LGKK_SCOPE_COARSE && mishLayer == 0)
        return false;
    return true;
}

void DrawGraph::clear_lgkk_level_cache()
{
    lgkkCacheActiveCount = 0;
    lgkkCacheMisfLevel = -1;
    lgkkCacheScaleL0 = 1.0;
    lgkkActiveIndex.clear();
    lgkkDistanceMatrix.clear();
    lgkkPairs.clear();
}

//****************************************************************
//
//    memory exception handler
//
//****************************************************************
void DrawGraph::noMoreMemory(){
        throw std::bad_alloc();
}

/**
 * @file DrawGraph.cpp
 * @brief Implementation of the Maximal Independent Set Filtration (MISF) algorithm.
 */

/**
 * @brief Constructs a Maximal Independent Set Filtration (MISF) for the graph.
 *
 * This function implements the MISF algorithm, which creates a hierarchical
 * decomposition of the graph into nested independent sets. Each level of the
 * filtration contains vertices that are increasingly far apart.
 *
 * @details The MISF algorithm works as follows:
 * 1. Start with the full vertex set V_0 = V.
 * 2. For each level i (starting from 1):
 *    a. Randomly select vertices from V_(i-1) to form V_i.
 *    b. Ensure that the graph distance between any two vertices in V_i is at least 2^(i-1) + 1.
 *    c. Continue until V_i is maximal (no more vertices can be added).
 * 3. Repeat steps 2a-c for increasing i until the filtration depth is reached.
 *
 * The implementation uses efficient data structures and algorithms to achieve
 * an average time complexity of O(n log n) for most graphs.
 *
 * @pre The graph structure (adjacency lists, etc.) must be initialized before calling this function.
 * @post The MISF structure is constructed and stored in the class member variables.
 *
 * @note This function modifies the graph structure in-place.
 *
 * @see bfs_cmisf() for the breadth-first search implementation used in distance calculations.
 *
 * @todo Parallelize BFS computations for improved performance on large graphs.
 * @todo Investigate more sophisticated random sampling techniques to improve independent set quality.
 *
 * Key member variables used/modified:
 * @li mish: Array representing the MISF structure.
 * @li misfSize: Array storing the size of each MISF level.
 * @li vertDepth: Array recording the depth (level) of each vertex in the filtration.
 * @li numOfVert: Total number of vertices in the graph.
 * @li numOfInitVert: Number of initial vertices to consider.
 * @li misfLevel: The current level of MISF being constructed.
 * @li initMishHeight: The final height of the MISF structure.
 *
 * @throw std::bad_alloc If memory allocation for bfsVectQueue fails.
 * @throw std::runtime_error If assertions on misfLevel or log_2_n fail.
 *
 * @warning This function assumes that the graph is connected. Behavior for disconnected graphs is undefined.
 */
void DrawGraph::create_misf() {
    // Initialize variables for MISF construction
    size_tt depthLim;
    size_tt shift = 0;
    std::queue<size_tt> **bfsVectQueue = new std::queue<size_tt>*[numOfVert];
    for(size_tt i = 0; i < numOfVert; i++)
        bfsVectQueue[i] = nullptr;

    // Start with the first level of MISF - a standard independent subset of vertices of the graph
    misfLevel = 1;
    misfSize[misfLevel-1] = numOfVert;
    size_tt mishSizeCurrLevel = numOfVert;

    size_tt itr;
    do {
        // Prepare for the next level of MISF
        size_tt mishSizePrevLevel = mishSizeCurrLevel;
        mishSizeCurrLevel = 0;
        itr = 0;     // independent set iterator

        // Construct the current level of MISF
        while (mishSizePrevLevel > numOfInitVert && itr < mishSizePrevLevel) {
            // Randomly select a vertex
            size_tt vert = itr + graph.fast_Rand() % (mishSizePrevLevel - itr);

            // Swap selected vertex to the current position
            std::swap(mish[vert], mish[itr]);
            inv[mish[vert]] = vert;
            inv[mish[itr]] = itr;

            // Add the vertex to the current MISF level
            std::swap(mish[mishSizeCurrLevel], mish[itr]);
            inv[mish[mishSizeCurrLevel]] = mishSizeCurrLevel;
            inv[mish[itr]] = itr;
            itr++;

            // Process the newly added vertex
            size_tt newEl = mish[mishSizeCurrLevel++];
            vertDepth[newEl] = misfLevel;

            // Remove vertices within the required distance
            if (misfLevel == 1) {
                // For the first level, remove direct neighbors
                for (size_tt i = 0; i < deg[newEl] && itr < numOfVert; i++) {
                    size_tt adj = graph.adjList[newEl+1][i];
                    if (inv[adj] >= itr) {
                        std::swap(mish[itr], mish[inv[adj]]);
                        inv[mish[inv[adj]]] = inv[adj];
                        inv[mish[itr]] = itr;
                        itr++;
                    }
                }
            } else {
                // For higher levels, use BFS results to remove vertices
                for (size_tt shiftedDepth = 0; shiftedDepth < depthLim-shift+1; shiftedDepth++) {
                    while (!bfsVectQueue[newEl][shiftedDepth].empty()) {
                        size_tt adj = bfsVectQueue[newEl][shiftedDepth].front();
                        bfsVectQueue[newEl][shiftedDepth].pop();
                        if (inv[adj] >= itr) {
                            std::swap(mish[itr], mish[inv[adj]]);
                            inv[mish[inv[adj]]] = inv[adj];
                            inv[mish[itr]] = itr;
                            itr++;
                        }
                    }
                }
            }
        }

        // Prepare for the next level
        depthLim = (size_tt)pow(2, misfLevel);
        depthLim = std::min(depthLim, (size_tt)(numOfVert-1));

        assert(misfLevel < log_2_n);
        misfSize[misfLevel] = mishSizeCurrLevel;

        // Prepare BFS queues for the next level
        for (int i = 0; i < numOfVert; i++)
            marked[i] = -1;
        shift = (size_tt)pow(2, misfLevel-1) + 1;
        for (size_tt i = 0; i < mishSizeCurrLevel; i++) {
            if(bfsVectQueue[mish[i]]){
                delete [] bfsVectQueue[mish[i]];
                bfsVectQueue[mish[i]] = nullptr;
            }
            bfsVectQueue[mish[i]] = new std::queue<size_tt>[depthLim-shift+1];
            bfs_cmisf(mish[i], bfsVectQueue[mish[i]], shift, depthLim);
        }
        misfLevel++;
    } while (itr);

    // Finalize MISF construction
    assert(misfLevel <= log_2_n);
    misfLevel -= 1;
    while (misfSize[misfLevel] < numOfInitVert)
        misfLevel--;

    size_tt v = 0;
    if (misfSize[misfLevel] > numOfInitVert) {
        while (v < numOfInitVert)
            vertDepth[mish[v++]] = misfLevel + 1;
        misfLevel++;
        misfSize[misfLevel] = numOfInitVert;
    }

    //debug("misfLevel=" << misfLevel);
    initMishHeight = misfLevel;

    // Clean up
    for (size_tt i = 0; i < numOfVert; i++){
        if(bfsVectQueue[i])
            delete [] bfsVectQueue[i];
    }
    delete [] bfsVectQueue;
}

/**
 * @brief Performs a modified Breadth-First Search for Maximal Independent Set Filtration.
 *
 * This function implements a specialized version of Breadth-First Search (BFS)
 * used in the creation of a Maximal Independent Set Filtration (MISF). It explores
 * the graph from a given root vertex up to a specified depth limit, collecting
 * vertices that are relevant to the current MISF level.
 *
 * @param root The starting vertex for the BFS.
 * @param bfsVectQueue An array of queues to store the BFS results for each depth level.
 * @param shift The offset used to adjust the depth at which vertices are added to bfsVectQueue.
 * @param depthLim The maximum depth to explore in the BFS.
 *
 * @pre The graph structure must be initialized.
 * @pre The MISF process must have begun (misfLevel and misfSize must be set).
 * @pre marked array must be initialized with values != root for all vertices.
 *
 * @post marked array will have all explored vertices marked with the root's index.
 * @post bfsVectQueue will contain vertices relevant to the MISF process at appropriate depths.
 *
 * @details The function performs the following steps:
 * 1. Initializes the BFS from the root vertex.
 * 2. Explores adjacent vertices up to the depth limit.
 * 3. Marks visited vertices to avoid revisiting.
 * 4. Adds relevant vertices (within the current MISF level size) to bfsVectQueue.
 * 5. Terminates when the depth limit is reached or no more vertices are available to explore.
 *
 * @note This function is specifically designed for use in the MISF algorithm and
 *       may not be suitable for general-purpose BFS applications.
 *
 * @see create_misf() for the main MISF algorithm that uses this function.
 *
 * @warning This function modifies the global 'marked' array. Ensure proper
 *          initialization before calling this function.
 *
 * Key member variables used:
 * @li deg: Array of vertex degrees.
 * @li graph.adjList: Adjacency list representation of the graph.
 * @li inv: Inverse mapping of vertex indices.
 * @li marked: Array to mark visited vertices.
 * @li misfSize: Array storing the size of each MISF level.
 * @li misfLevel: The current level of MISF being constructed.
 * @li numOfVert: Total number of vertices in the graph.
 *
 * @throw std::bad_alloc If memory allocation for FastQueue fails.
 *
 * Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges
 * within the depth limit from the root.
 * Space Complexity: O(V) for the FastQueue and additional data structures.
 */
void DrawGraph::bfs_cmisf(size_tt root,
                          std::queue<size_tt> *bfsVectQueue,
                          size_tt shift,
                          size_tt depthLim)
{
    size_tt vert = root;
    size_tt currentDepth = 1;
    FastQueue< size_tt > vertDepthQueue(4*numOfVert); // an array based queue
    // consisitaing of vertex, its dist from root, vertex, its dist from root,...

    marked[root] = root; // all elements discovered in this BFS are marked
    // with index root

    size_tt mishSizeLim = misfSize[misfLevel];

    while( currentDepth <= depthLim ){
        // Mark and process vert's un_marked adjacent vertices,
        // and place them in the queue.
        for (size_tt adjVert = 0; adjVert < deg[vert]; adjVert++){
            size_tt overt = graph.adjList[vert+1][adjVert];
            if ( marked[overt] != (int)root){
                marked[overt] = root;
                vertDepthQueue.enqueue(overt);
                vertDepthQueue.enqueue(currentDepth);

                if(inv[overt] < mishSizeLim)
                    bfsVectQueue[currentDepth-shift].push(overt);
            }
        }
        if(!vertDepthQueue.is_empty( )){
            vert = vertDepthQueue.dequeue();
            currentDepth = vertDepthQueue.dequeue() + 1;
        } else
            break;
    }
}
