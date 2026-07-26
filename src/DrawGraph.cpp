#include "DrawGraph.h"
#include <algorithm>
#include <functional>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <limits>
#include <sstream>

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
                     size_tt _lgkkRoundsCoarse,
                     size_tt _lgkkRoundsPreFinal,
                     size_tt _lgkkRoundsFinal,
                     size_tt _lgkkLocalNbrs,
                     size_tt _lgkkLandmarkCount,
                     size_tt _lgkkScope,
                     size_tt _lgkkActiveLimit,
                     bool _weightedCore)
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
  metricNeighborCap(0),
  metricScratchEpoch(0),
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
  lgkkRoundsCoarse(_lgkkRoundsCoarse),
  lgkkRoundsPreFinal(_lgkkRoundsPreFinal),
  lgkkRoundsFinal(_lgkkRoundsFinal),
  lgkkLocalNbrs(_lgkkLocalNbrs),
  lgkkLandmarkCount(_lgkkLandmarkCount),
  lgkkScope((_lgkkScope == LGKK_SCOPE_COARSE)
                ? LGKK_SCOPE_COARSE
                : LGKK_SCOPE_ALL),
  lgkkActiveLimit(std::max<size_tt>(1, _lgkkActiveLimit)),
  weightedCore(_weightedCore),
  activeVertCount(0),
  misfLevel(0),
  initMishHeight(0),
  currentRoundInLevel(0),
  finalAnchorReady(false),
  traceMode(TRACE_NONE),
  traceEvery(1),
  traceLevelIndex(0),
  finalAnchorPos(_graph.get_numOfVert()),
  refinementStepTraceEnabled(false),
  collectRefinementStepDetails(false),
  refinementStepTraceLevelIndex(-1),
  refinementStepTraceMisfLevel(-1),
  refinementStepTraceRoundStart(-1),
  refinementStepTraceRoundEnd(-1),
  insertionTraceEnabled(false),
  lastAttractionDisp(static_cast<size_t>(_dim), 0.0),
  lastRepulsionDisp(static_cast<size_t>(_dim), 0.0),
  lgkkCacheActiveCount(0),
  lgkkCacheMisfLevel(-1),
  lgkkCacheScaleL0(1.0)
{
    if( dim != 2 && dim != 3 )
        throw std::runtime_error("only 2D and 3D layouts are supported");

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
    for(size_tt i = 0; i < numOfVert; i++)
        marked[i] = -1;
    
    log_2_n = ilog(numOfVert)+2;
    misfSize  = new size_tt[log_2_n];
    vertDepth = new size_tt[numOfVert];
    for(size_tt i = 0; i < log_2_n; i++)
        misfSize[i] = 0;
    for(size_tt i = 0; i < numOfVert; i++)
        vertDepth[i] = 0;

    if(numOfVert == numOfInitVert){
        misfSize[0] = numOfVert;
        misfLevel = 0;
    } else if(weightedCore)
        create_misf_weighted();
    else
        create_misf();
    initMishHeight = misfLevel;
    
    // resetting marked array values for BFSs in the main part of
    // the program - MishEngine(s).cpp
    for(size_tt i = 0; i < numOfVert; i++){
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
        

    // numOfNbrs directly controls the retained local graph-distance
    // neighborhood size used by the refinement routines.
    nbr       = new size_tt[log_2_n];
    itr = 0;
    while( itr < log_2_n && misfSize[itr] ){
        size_tt available = misfSize[itr] > 0 ? misfSize[itr] - 1 : 0;
        nbr[itr] = std::min(numOfNbrs, available);
        itr++;
    }

    nbrs = new size_tt**[numOfVert];
    nbrsDepth = new size_tt[numOfVert];
    metricNbrs = new std::vector<MetricNeighbor>*[numOfVert];
    metricNbrsDepth = new size_tt[numOfVert];
    for(size_tt i = 0; i < numOfVert; i++){
        nbrs[i] = nullptr;
        nbrsDepth[i] = 0;
        metricNbrs[i] = nullptr;
        metricNbrsDepth[i] = 0;
    }
    diam = 0;
    double metricHeight = 0.0;
    size_tt height = 0;
    for(size_tt i = 0; i < numOfInitVert; i++){
        if(weightedCore){
            metricHeight = metric_me_init_v1(mish[i]);
            size_tt roundedHeight = static_cast<size_tt>(std::ceil(metricHeight));
            if( diam < roundedHeight )
                diam = roundedHeight;
        } else {
            height = bfs_me_init_v2(mish[i]);
            if( diam < height )
                diam = height;
        }
    }
        
    boxSize = (coord_t)(edge * .7 * diam);
    box2Size = 2 * boxSize + 1;
        
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

void DrawGraph::configure_refinement_step_trace(int level_index,
                                                int misf_level,
                                                int round_start,
                                                int round_end)
{
    refinementStepTraceEnabled = true;
    refinementStepTraceLevelIndex = level_index;
    refinementStepTraceMisfLevel = misf_level;
    refinementStepTraceRoundStart = round_start;
    refinementStepTraceRoundEnd = std::max(round_start, round_end);
    refinementStepTrace = RefinementStepTrace();
}

void DrawGraph::configure_insertion_trace(bool enabled)
{
    insertionTraceEnabled = enabled;
    insertionTrace = InsertionTrace();
}

void DrawGraph::configure_weighted_metric_search(size_tt maxSettled)
{
    metricNeighborCap = maxSettled;
}

bool DrawGraph::should_record_refinement_step(size_tt level_index,
                                              size_tt misf_level,
                                              size_tt round) const
{
    return refinementStepTraceEnabled &&
        static_cast<int>(level_index) == refinementStepTraceLevelIndex &&
        static_cast<int>(misf_level) == refinementStepTraceMisfLevel &&
        static_cast<int>(round) >= refinementStepTraceRoundStart &&
        static_cast<int>(round) <= refinementStepTraceRoundEnd;
}

size_t DrawGraph::record_refinement_step_pre(
    size_tt vert,
    size_tt order_index,
    size_tt active_count,
    size_tt level_index,
    size_tt misf_level,
    size_tt round,
    const Point<> &coord_before,
    const Point<> &pre_temp_disp,
    const std::vector<double> &attraction_disp,
    const std::vector<double> &repulsion_disp,
    const Point<> &applied_disp,
    coord_t heat_before,
    coord_t heat_after,
    float old_cos_before,
    float old_cos_after,
    coord_t old_disp_norm_before,
    coord_t pre_temp_disp_norm,
    const std::string &attraction_edges,
    const std::string &repulsion_neighbors)
{
    std::vector<double> coord_before_flat(dim);
    std::vector<double> pre_temp_disp_flat(dim);
    std::vector<double> attraction_disp_flat(dim);
    std::vector<double> repulsion_disp_flat(dim);
    std::vector<double> applied_disp_flat(dim);
    std::vector<double> coord_after_flat(dim,
                                         std::numeric_limits<double>::quiet_NaN());
    coord_before_flat[0] = coord_before.getX();
    pre_temp_disp_flat[0] = pre_temp_disp.getX();
    attraction_disp_flat[0] = attraction_disp[0];
    repulsion_disp_flat[0] = repulsion_disp[0];
    applied_disp_flat[0] = applied_disp.getX();
    if(dim > 1){
        coord_before_flat[1] = coord_before.getY();
        pre_temp_disp_flat[1] = pre_temp_disp.getY();
        attraction_disp_flat[1] = attraction_disp[1];
        repulsion_disp_flat[1] = repulsion_disp[1];
        applied_disp_flat[1] = applied_disp.getY();
    }
    if(dim > 2){
        coord_before_flat[2] = coord_before.getZ();
        pre_temp_disp_flat[2] = pre_temp_disp.getZ();
        attraction_disp_flat[2] = attraction_disp[2];
        repulsion_disp_flat[2] = repulsion_disp[2];
        applied_disp_flat[2] = applied_disp.getZ();
    }

    refinementStepTrace.level_indices.push_back(static_cast<int>(level_index));
    refinementStepTrace.misf_levels.push_back(static_cast<int>(misf_level));
    refinementStepTrace.rounds.push_back(static_cast<int>(round));
    refinementStepTrace.active_counts.push_back(static_cast<int>(active_count));
    refinementStepTrace.order_indices.push_back(static_cast<int>(order_index));
    refinementStepTrace.vertices.push_back(static_cast<int>(vert) + 1);
    refinementStepTrace.heat_before.push_back(heat_before);
    refinementStepTrace.heat_after.push_back(heat_after);
    refinementStepTrace.old_cos_before.push_back(old_cos_before);
    refinementStepTrace.old_cos_after.push_back(old_cos_after);
    refinementStepTrace.old_disp_norm_before.push_back(old_disp_norm_before);
    refinementStepTrace.pre_temp_disp_norm.push_back(pre_temp_disp_norm);
    refinementStepTrace.coords_before.push_back(std::move(coord_before_flat));
    refinementStepTrace.pre_temp_disp.push_back(std::move(pre_temp_disp_flat));
    refinementStepTrace.attraction_disp.push_back(std::move(attraction_disp_flat));
    refinementStepTrace.repulsion_disp.push_back(std::move(repulsion_disp_flat));
    refinementStepTrace.applied_disp.push_back(std::move(applied_disp_flat));
    refinementStepTrace.coords_after.push_back(std::move(coord_after_flat));
    refinementStepTrace.attraction_edges.push_back(attraction_edges);
    refinementStepTrace.repulsion_neighbors.push_back(repulsion_neighbors);
    const int parentRow = static_cast<int>(refinementStepTrace.vertices.size());
    for(size_t i = 0; i < lastAttractionTermNeighbors.size(); i++){
        refinementStepTrace.attraction_term_parent_rows.push_back(parentRow);
        refinementStepTrace.attraction_term_indices.push_back(static_cast<int>(i) + 1);
        refinementStepTrace.attraction_term_vertices.push_back(static_cast<int>(vert) + 1);
        refinementStepTrace.attraction_term_neighbors.push_back(lastAttractionTermNeighbors[i]);
        refinementStepTrace.attraction_term_weights.push_back(lastAttractionTermWeights[i]);
        refinementStepTrace.attraction_term_norm2.push_back(lastAttractionTermNorm2[i]);
        refinementStepTrace.attraction_term_desired.push_back(lastAttractionTermDesired[i]);
        refinementStepTrace.attraction_term_desired2.push_back(lastAttractionTermDesired2[i]);
        refinementStepTrace.attraction_term_scale.push_back(lastAttractionTermScale[i]);
        refinementStepTrace.attraction_term_delta.push_back(lastAttractionTermDelta[i]);
        refinementStepTrace.attraction_term_step.push_back(lastAttractionTermStep[i]);
        refinementStepTrace.attraction_term_cumulative.push_back(lastAttractionTermCumulative[i]);
    }
    return refinementStepTrace.vertices.size() - 1;
}

void DrawGraph::record_refinement_step_after(size_t row, size_tt vert)
{
    if(row >= refinementStepTrace.coords_after.size())
        return;
    refinementStepTrace.coords_after[row][0] = pos[vert].getX();
    if(dim > 1)
        refinementStepTrace.coords_after[row][1] = pos[vert].getY();
    if(dim > 2)
        refinementStepTrace.coords_after[row][2] = pos[vert].getZ();
}

void DrawGraph::record_insertion_trace(
    size_tt root,
    size_tt level_index,
    size_tt misf_level_arg,
    size_tt previous_active_count,
    size_tt active_count,
    size_tt root_depth,
    size_tt anchor_count_requested,
    size_tt anchor_count_used,
    size_tt insertion_mode,
    size_tt local_kk_steps,
    const std::vector<size_tt> &anchors,
    const std::vector<double> &anchor_dist,
    const Point<> &coord_initial,
    const Point<> &coord_after,
    const Point<> &old_disp_initial,
    const Point<> &old_disp_after,
    coord_t old_disp_norm_initial_arg,
    coord_t old_disp_norm_after_arg)
{
    if(!insertionTraceEnabled)
        return;

    auto flatten_point = [&](const Point<> &point) {
        std::vector<double> out(dim, 0.0);
        out[0] = point.getX();
        if(dim > 1)
            out[1] = point.getY();
        if(dim > 2)
            out[2] = point.getZ();
        return out;
    };

    std::ostringstream anchor_stream;
    anchor_stream << std::setprecision(17);
    const size_tt limit = std::min<size_tt>(
        std::min<size_tt>(anchor_count_used, anchors.size()),
        anchor_dist.size());
    for(size_tt i = 0; i < limit; i++){
        if(i > 0)
            anchor_stream << ";";
        anchor_stream << static_cast<int>(anchors[i]) + 1 << ":"
                      << anchor_dist[i];
    }

    insertionTrace.level_indices.push_back(static_cast<int>(level_index));
    insertionTrace.misf_levels.push_back(static_cast<int>(misf_level_arg));
    insertionTrace.previous_active_counts.push_back(
        static_cast<int>(previous_active_count));
    insertionTrace.active_counts.push_back(static_cast<int>(active_count));
    insertionTrace.order_indices.push_back(static_cast<int>(inv[root]) + 1);
    insertionTrace.vertices.push_back(static_cast<int>(root) + 1);
    insertionTrace.root_depths.push_back(static_cast<int>(root_depth));
    insertionTrace.anchor_count_requested.push_back(
        static_cast<int>(anchor_count_requested));
    insertionTrace.anchor_count_used.push_back(static_cast<int>(anchor_count_used));
    insertionTrace.insertion_modes.push_back(static_cast<int>(insertion_mode));
    insertionTrace.local_kk_steps.push_back(static_cast<int>(local_kk_steps));
    insertionTrace.anchors.push_back(anchor_stream.str());
    insertionTrace.old_disp_norm_initial.push_back(old_disp_norm_initial_arg);
    insertionTrace.old_disp_norm_after.push_back(old_disp_norm_after_arg);
    insertionTrace.coords_initial.push_back(flatten_point(coord_initial));
    insertionTrace.coords_after.push_back(flatten_point(coord_after));
    insertionTrace.old_disp_initial.push_back(flatten_point(old_disp_initial));
    insertionTrace.old_disp_after.push_back(flatten_point(old_disp_after));
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
    if(metricNbrs){
        for(size_tt i = 0; i < numOfVert; i++){
            if(metricNbrs[i])
                delete [] metricNbrs[i];
        }
        delete [] metricNbrs;
    }
    delete [] metricNbrsDepth;

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
        std::vector<size_tt> selectedDist;
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

size_tt DrawGraph::lgkk_round_budget_for_layer(size_tt mishLayer) const
{
    size_tt budget = 0;
    if(mishLayer == 0)
        budget = lgkkRoundsFinal;
    else if(mishLayer == 1)
        budget = lgkkRoundsPreFinal;
    else
        budget = lgkkRoundsCoarse;

    if(budget == 0)
        budget = lgkkMultiscaleRounds;
    return budget;
}

bool DrawGraph::should_run_multiscale_lgkk(size_tt activeCount,
                                           size_tt mishLayer) const
{
    if(lgkk_round_budget_for_layer(mishLayer) == 0)
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

/**
 * Construct a maximal independent set filtration for the graph.
 *
 * Each successive level is a maximal subset whose vertices satisfy the
 * level-specific graph-distance separation.
 */
void DrawGraph::create_misf() {
    // Initialize variables for MISF construction
    size_tt depthLim = 0;
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

        if(misfLevel >= log_2_n)
            throw std::runtime_error("MISF level exceeded the allocated hierarchy");
        misfSize[misfLevel] = mishSizeCurrLevel;

        // Prepare BFS queues for the next level
        for (size_tt i = 0; i < numOfVert; i++)
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
    if(misfLevel > log_2_n)
        throw std::runtime_error("MISF level exceeded the allocated hierarchy");
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
    // consisting of vertex, its distance from root, vertex, its distance from root, ...

    marked[root] = root; // all elements discovered in this BFS are marked
    // with index root

    size_tt mishSizeLim = misfSize[misfLevel];

    while( currentDepth <= depthLim ){
        // Mark and process vert's unmarked adjacent vertices.
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
