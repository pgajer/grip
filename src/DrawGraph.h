// DrawGraph.hpp a header file for DrawGraph class

#ifndef DRAW_GRAPH_HPP
#define DRAW_GRAPH_HPP

#include <iostream>
#include <new>
#include <cmath>
#include <queue>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include <functional>

#include "Point.h"
#include "Graph.h"
#include "FastQueue.h"
#include "Debug.h"

using size_tt = uint32_t;

#define ROUND(a)       ((a)>0 ? (int)((a)+0.5) : -(int)(0.5-(a)))
#define ROUND_L(a) ((a)>0 ? (unsigned long)((a)+0.5) : -(unsigned long)(0.5-(a)))
#define PLACEMENT_BARYCENTER 0
#define PLACEMENT_CIRCLE 1
#define TRACE_NONE 0
#define TRACE_ROUND 1
#define TRACE_LEVEL 2
#define FINAL_STAGE_FR 0
#define FINAL_STAGE_KK_REPULSE 1
#define LEVEL0_INSERT_INHERIT 100
#define LEVEL0_INSERT_BARYCENTER 101
#define LEVEL0_INSERT_LEAST_SQUARES 102
#define INSERT_ANCHOR_SCOPE_ANY_HIGHER 0
#define INSERT_ANCHOR_SCOPE_PREV_MISF 1
#define INSERT_ANCHOR_STRATEGY_FIRST 0
#define INSERT_ANCHOR_STRATEGY_DISTANCE_BAND 1
#define INSERT_ANCHOR_STRATEGY_BALANCED_BAND 2
#define INSERT_ANCHOR_STRATEGY_SPREAD_PREV 3
#define LGKK_SCOPE_ALL 0
#define LGKK_SCOPE_COARSE 1

//**************************************************************
//
//	Class name : DrawGraph - compact version
//      for the GRIP layout engine
//
//**************************************************************
class DrawGraph 
{
  public:
    struct MetricNeighbor {
        size_tt vert;
        double dist;
    };
    struct RefinementStepTrace {
        std::vector<int> level_indices;
        std::vector<int> misf_levels;
        std::vector<int> rounds;
        std::vector<int> active_counts;
        std::vector<int> order_indices;
        std::vector<int> vertices;
        std::vector<double> heat_before;
        std::vector<double> heat_after;
        std::vector<double> old_cos_before;
        std::vector<double> old_cos_after;
        std::vector<double> old_disp_norm_before;
        std::vector<double> pre_temp_disp_norm;
        std::vector<std::vector<double>> coords_before;
        std::vector<std::vector<double>> pre_temp_disp;
        std::vector<std::vector<double>> attraction_disp;
        std::vector<std::vector<double>> repulsion_disp;
        std::vector<std::vector<double>> applied_disp;
        std::vector<std::vector<double>> coords_after;
        std::vector<std::string> attraction_edges;
        std::vector<std::string> repulsion_neighbors;
        std::vector<int> attraction_term_parent_rows;
        std::vector<int> attraction_term_indices;
        std::vector<int> attraction_term_vertices;
        std::vector<int> attraction_term_neighbors;
        std::vector<double> attraction_term_weights;
        std::vector<double> attraction_term_norm2;
        std::vector<double> attraction_term_desired;
        std::vector<double> attraction_term_desired2;
        std::vector<double> attraction_term_scale;
        std::vector<std::vector<double>> attraction_term_delta;
        std::vector<std::vector<double>> attraction_term_step;
        std::vector<std::vector<double>> attraction_term_cumulative;
    };
    struct InsertionTrace {
        std::vector<int> level_indices;
        std::vector<int> misf_levels;
        std::vector<int> previous_active_counts;
        std::vector<int> active_counts;
        std::vector<int> order_indices;
        std::vector<int> vertices;
        std::vector<int> root_depths;
        std::vector<int> anchor_count_requested;
        std::vector<int> anchor_count_used;
        std::vector<int> insertion_modes;
        std::vector<int> local_kk_steps;
        std::vector<std::string> anchors;
        std::vector<double> old_disp_norm_initial;
        std::vector<double> old_disp_norm_after;
        std::vector<std::vector<double>> coords_initial;
        std::vector<std::vector<double>> coords_after;
        std::vector<std::vector<double>> old_disp_initial;
        std::vector<std::vector<double>> old_disp_after;
    };

      friend class MesaPlot;
    
    DrawGraph(const Graph &_graph,
              size_tt _dim         ,
              size_tt _rounds      ,
              size_tt _finalRounds ,
              size_tt _tinit_factor,
              size_tt _numOfInitVert,
              size_tt _numOfNbrs   ,
              coord_t _r             ,  //parameters of update_Local_Temp_v3()
              coord_t _s             ,
              coord_t _repulsionFactor,
              size_tt _placementMode,
              bool _displayPar     ,
              size_tt _finalStageMode = FINAL_STAGE_FR,
              coord_t _coarseRepulsionFactor = 0.0,
              size_tt _coarseRepulsionSample = 0,
              size_tt _coarseRepulsionExactBelow = 0,
              coord_t _finalAnchorFactor = 0.0,
              coord_t _finalMoveScaleAfterFirst = 1.0,
              size_tt _insertionAnchorCount = 3,
              size_tt _insertionAnchorScope = INSERT_ANCHOR_SCOPE_ANY_HIGHER,
              size_tt _insertionAnchorStrategy = INSERT_ANCHOR_STRATEGY_FIRST,
              size_tt _level0InsertionMode = LEVEL0_INSERT_INHERIT,
              size_tt _level0AnchorCount = 3,
              size_tt _level0LocalKkSteps = 3,
              size_tt _lgkkMultiscaleRounds = 0,
              size_tt _lgkkRoundsCoarse = 0,
              size_tt _lgkkRoundsPreFinal = 0,
              size_tt _lgkkRoundsFinal = 0,
              size_tt _lgkkLocalNbrs = 20,
              size_tt _lgkkLandmarkCount = 8,
              size_tt _lgkkScope = LGKK_SCOPE_ALL,
              size_tt _lgkkActiveLimit = 4096,
              bool _weightedCore = false);
    
    ~DrawGraph();

    //
    // GRIP layout core
    //
    size_tt bfs_me_init_v2(size_tt root);
    void bfs_me_v4(size_tt root);
    double metric_me_init_v1(size_tt root);
    void metric_me_v1(size_tt root);
    void bfs_cmisf(size_tt root,
                   std::queue<size_tt> *bfsVectQueue,
                   size_tt shift,   // bfsVector starts at dist shift from root
                   size_tt depthLim);
    void create_misf();  // version of create_mish() creating also vertDepth
    void create_misf_weighted();
    void KK_spring_v4(const size_tt root,
                      size_tt *rootNbrsLayer,
                      size_tt mishLayer);
    void KK_spring_weighted_v1(const size_tt root,
                               const std::vector<MetricNeighbor> &rootNbrsLayer,
                               size_tt mishLayer);
    void KK_spring_final(const size_tt root,
                         size_tt *rootNbrsLayer,
                         size_tt mishLayer);
    void KK_spring_weighted_final_v1(const size_tt root,
                                     const std::vector<MetricNeighbor> &rootNbrsLayer,
                                     size_tt mishLayer);
    void FR_spring_v2(const size_tt vert,
                      size_tt *vertNbrs,
                      size_tt misfLayer);
    void FR_spring_weighted_v1(const size_tt vert,
                               const std::vector<MetricNeighbor> &vertNbrs,
                               size_tt misfLayer);
    void KK_spring_local(const size_tt vert,
                         size_tt *closeVert,
                         size_tt *closeVertDist,
                         size_tt size);
    void KK_spring_weighted_local_v1(const size_tt vert,
                                     const size_tt *closeVert,
                                     const double *closeVertDist,
                                     size_tt size);
    void update_Local_Temp_v2( size_tt vert );
    void update_Local_Temp_v3( size_tt vert, coord_t r, coord_t s );
    
    void mish_engine();
    void mish_engine_weighted();
    // FR_spring() is used in the GRIP layout core and utilizes
    // attractive/repulsive force schedule
    void FR_spring(const size_tt root,
                   size_tt *rootNbrsLayer,
                   size_tt mishLayer);
    void add_coarse_global_repulsion(const size_tt vert,
                                     size_tt activeCount);
    void add_active_global_repulsion(const size_tt vert,
                                     size_tt activeCount,
                                     coord_t repulsionScale);
    void add_active_global_repulsion_exact(const size_tt vert,
                                           size_tt activeCount,
                                           coord_t repulsionScale);
    void add_active_global_repulsion_sampled(const size_tt vert,
                                             size_tt activeCount,
                                             size_tt sampleCount,
                                             coord_t repulsionScale);
    void add_coarse_global_repulsion_exact(const size_tt vert,
                                           size_tt activeCount);
    void add_coarse_global_repulsion_sampled(const size_tt vert,
                                             size_tt activeCount,
                                             size_tt sampleCount);
    void add_final_anchor_force(const size_tt vert);
    

    // Generate a list of random positions of numOfVert points in dim
    // dimensional centered at 0 cube whose edge has length = width
    // parameter
    void rand_Positions();
    void print_Positions(std::ostream &output = std::cout);
    void print_Positions_To_File(const char *posFile = NULL);
    void read_Positions_From_File(const char *posFile);

    size_tt get_Dim() const { return dim; }
    size_tt get_Diam() const { return diam; }
    size_tt get_NumOfVert() const { return numOfVert; }
    Point<> *get_Pos() const { return pos; }
    void configure_trace(size_tt mode, size_tt every);
    void configure_refinement_step_trace(int level_index,
                                         int misf_level,
                                         int round_start,
                                         int round_end);
    void configure_insertion_trace(bool enabled);
    void configure_weighted_metric_search(size_tt maxSettled);
    const std::vector<std::vector<double>> &get_trace_frames() const { return traceFrames; }
    const std::vector<std::string> &get_trace_phases() const { return tracePhases; }
    const std::vector<int> &get_trace_level_indices() const { return traceLevelIndices; }
    const std::vector<int> &get_trace_misf_levels() const { return traceMisfLevels; }
    const std::vector<int> &get_trace_rounds() const { return traceRounds; }
    const std::vector<int> &get_trace_active_counts() const { return traceActiveCounts; }
    const RefinementStepTrace &get_refinement_step_trace() const
    {
        return refinementStepTrace;
    }
    const InsertionTrace &get_insertion_trace() const
    {
        return insertionTrace;
    }
    size_tt get_MisfHeight() const { return initMishHeight; }
    size_tt get_MisfLevel() const { return misfLevel; }
    size_tt get_Mish(size_tt idx) const { return mish[idx]; }
    size_tt get_MisfSize(size_tt level) const { return misfSize[level]; }
    size_tt get_VertDepth(size_tt vert) const { return vertDepth[vert]; }
    size_tt get_NbrCount(size_tt level) const { return nbr[level]; }
    coord_t get_Edge(){ return edge; }
    coord_t get_Edge2(){ return edge2; }
    size_tt get_inv(size_tt vert){return inv[vert];}
    size_tt get_prevMishSize(){return prevMishSize;}
    coord_t dist(const Point<> &p, const Point<> &q);
    bool should_record_refinement_step(size_tt level_index,
                                       size_tt misf_level,
                                       size_tt round) const;
    size_t record_refinement_step_pre(size_tt vert,
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
                                      const std::string &repulsion_neighbors);
    void record_refinement_step_after(size_t row, size_tt vert);
    void record_insertion_trace(size_tt root,
                                size_tt level_index,
                                size_tt misf_level,
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
                                coord_t old_disp_norm_initial,
                                coord_t old_disp_norm_after);
    
    //memory exception handler
    static void noMoreMemory();

    bool createList;
    
private:
    const Graph &graph;               // reference to the graph we want to draw
    size_tt dim;                       // dimension of the drawing
    size_tt rounds;                       // max number of rounds [4*|V|]
    size_tt initRounds;
    size_tt finalRounds; 
    size_tt numOfVert;                 // this is inherited from graph
    // we make it a data of DrawGraph to avoid frequent referals
    // to graph.get_numOfVert()

    const coord_t edge;               // desired edge length [128]
    const coord_t edge2;              // = edge * edge
    coord_t tinit;              // desired minimal temp [3]
    const coord_t SMALL_DIST;         // abs value of disturbance [32]
    const coord_t SMALL_DIST2;        // = 2*SMALL_DIST + 1
    coord_t boxSize;   // width of the view window is 2*boxSize
    coord_t box2Size;                 // 2*boxSize + 1

    Point<> center;                   // = \sum pos[vert]
    Point<> baricenter;               // center / numOfVert
    Point<> *pos;        // positions of vertices vector
    Point<> *disp;       // displacement vector
    Point<> *oldDisp;    // old displacement vector
    coord_t *dispNorm;
    coord_t *oldDispNorm;
    int *marked;             // array of discovered by BFS vertices
    

    size_tt *deg;             // array of degrees
    float *old_cos;
    int *maxNorm;              // vector of maximal norms of disp
    coord_t *heat;                   // local temperature

    //
    // variable used specifically in mish_engines
    //
    size_tt *mish;               // maximal independent set hierarchy
    size_tt *inv;                // "inverse" of mish
    size_tt *misfSize;  // filtration sizes used by the GRIP layout core
    size_tt log_2_n;    // log_2(n)
    size_tt *vertDepth; // array of depths of vertices
    float AvgDeg;
    unsigned long maxCxty;
    unsigned long initCxty;
    size_tt smallLevel;
    size_tt *nbr;       // array of num of nbrs for each level of misf
    size_tt ***nbrs;    // storage of nbrs for corresponding levels
    size_tt *nbrsDepth; // number of allocated levels per vertex
    std::vector<MetricNeighbor> **metricNbrs; // weighted neighborhood caches
    size_tt *metricNbrsDepth;
    size_tt metricNeighborCap;
    std::vector<double> metricScratchDist;
    std::vector<uint32_t> metricScratchStamp;
    uint32_t metricScratchEpoch;
    size_tt prevSize;
    
    size_tt prevMishSize;
    size_tt diam;                         // graph's diameter
    size_tt numOfInitVert;
    size_tt numOfNbrs;
    coord_t r, s; // parameters of the main local temperature schedule
    bool listSwitch;
    bool displayPar;
    size_tt placementMode;
    size_tt finalStageMode;
    coord_t fedge2;
    coord_t coarseFedge2;
    size_tt coarseRepulsionSample;
    size_tt coarseRepulsionExactBelow;
    coord_t finalAnchorFactor;
    coord_t finalMoveScaleAfterFirst;
    size_tt insertionAnchorCount;
    size_tt insertionAnchorScope;
    size_tt insertionAnchorStrategy;
    size_tt level0InsertionMode;
    size_tt level0AnchorCount;
    size_tt level0LocalKkSteps;
    size_tt lgkkMultiscaleRounds;
    size_tt lgkkRoundsCoarse;
    size_tt lgkkRoundsPreFinal;
    size_tt lgkkRoundsFinal;
    size_tt lgkkLocalNbrs;
    size_tt lgkkLandmarkCount;
    size_tt lgkkScope;
    size_tt lgkkActiveLimit;
    bool weightedCore;
    size_tt activeVertCount;
    size_tt misfLevel;
    size_tt initMishHeight;
    size_tt currentRoundInLevel;
    bool finalAnchorReady;
    Point<> vect;    // an auxiliary vector holding vert1 - vert2
    size_tt traceMode;
    size_tt traceEvery;
    size_tt traceLevelIndex;
    std::vector<Point<>> finalAnchorPos;
    std::vector<std::vector<double>> traceFrames;
    std::vector<std::string> tracePhases;
    std::vector<int> traceLevelIndices;
    std::vector<int> traceMisfLevels;
    std::vector<int> traceRounds;
    std::vector<int> traceActiveCounts;
    bool refinementStepTraceEnabled;
    bool collectRefinementStepDetails;
    int refinementStepTraceLevelIndex;
    int refinementStepTraceMisfLevel;
    int refinementStepTraceRoundStart;
    int refinementStepTraceRoundEnd;
    RefinementStepTrace refinementStepTrace;
    bool insertionTraceEnabled;
    InsertionTrace insertionTrace;
    std::vector<double> lastAttractionDisp;
    std::vector<double> lastRepulsionDisp;
    std::string lastAttractionEdges;
    std::string lastRepulsionNeighbors;
    std::vector<int> lastAttractionTermNeighbors;
    std::vector<double> lastAttractionTermWeights;
    std::vector<double> lastAttractionTermNorm2;
    std::vector<double> lastAttractionTermDesired;
    std::vector<double> lastAttractionTermDesired2;
    std::vector<double> lastAttractionTermScale;
    std::vector<std::vector<double>> lastAttractionTermDelta;
    std::vector<std::vector<double>> lastAttractionTermStep;
    std::vector<std::vector<double>> lastAttractionTermCumulative;
    struct LgkkPathEdge {
        size_tt u;
        size_tt v;
    };
    struct LgkkPairCache {
        size_tt source;
        size_tt target;
        coord_t graphDistance;
        std::vector<LgkkPathEdge> pathEdges;
    };
    size_tt lgkkCacheActiveCount;
    int lgkkCacheMisfLevel;
    coord_t lgkkCacheScaleL0;
    std::vector<int> lgkkActiveIndex;
    std::vector<double> lgkkDistanceMatrix;
    std::vector<LgkkPairCache> lgkkPairs;

    //
    // SUPPORTING FUNCTIONS
    //

    // sched(x, ...) is constant for x <= max and x >= min and it is a linear
    // function connecting points ( max, maxVal) with (min, minVal) in the
    // interval (max, min)
    float sched(size_tt x,
                  size_tt max, size_tt maxVal, size_tt min, size_tt minVal){
        if( x <= max )
            return (float)maxVal;
        else if( max <= x && x <= min ){

            return ((minVal - maxVal)/(float)(min-max))*(x - max) + maxVal;
        } else
            return (float)minVal;

    }
        
    // schedule function with a parabolic piece
    size_tt sched2(size_tt x,
                   size_tt max, size_tt maxVal, size_tt min, size_tt minVal){
        if( x <= max )
            return maxVal;
        else if( max <= x && x <= min ){
            size_tt a = (size_tt)((minVal - maxVal)/(float)((max-min)*(max-min)));
            return a * (x - min)*(x - min) + minVal;
        } else
            return minVal;

    }

    // exponentially decaying schedule function
    // it assumes that we are working in the interval [0, min]
    // that is max = 0
    size_tt sched3(size_tt x,
                   size_tt max, size_tt maxVal, size_tt min, size_tt minVal){
        if( x <= max )
            return maxVal;
        else if( max <= x && x <= min ){
            double k = -log((double)minVal/maxVal)/min;
            
            return (size_tt)(ceil(maxVal * exp( -k * x )));
        } else
            return minVal;
    }

    double norm2(const Point<> &pt1, const Point<> &pt2){
        return ((pt1.getX() - pt2.getX()) * (pt1.getX() - pt2.getX()) +
                (pt1.getY() - pt2.getY()) * (pt1.getY() - pt2.getY()) +
                (pt1.getZ() - pt2.getZ()) * (pt1.getZ() - pt2.getZ()));
    }

    double norm2(size_tt u, size_tt v){
        return norm2(pos[u], pos[v]);
    }
    
    coord_t get_BoxSize(){ return boxSize; }

    coord_t sign( coord_t x )
        { 
        if( x > 0 )
            return 1.0;
        else if( x< 0 )
            return -1.0;
        else
            return 0.0;
        }

    coord_t sign( float x )
        { 
        if( x > 0 )
            return 1.0;
        else if( x < 0 )
            return -1.0;
        else
            return 0.0;
        }
    
    void update_Baricenter(){ baricenter = center /(coord_t)numOfVert; }
    Point<> get_Baricenter()
        {
            Point<> center;
            for(size_tt vert = 0; vert < numOfVert; vert++)
                center += pos[vert];
            return center /(coord_t)numOfVert;
        }
    // random point
    Point<> rand_Point();

    Point<> initial_position(const size_tt *closeVert,
                             const size_tt *closeVertDist,
                             size_tt count);
    Point<> initial_position_mode(const size_tt *closeVert,
                                  const size_tt *closeVertDist,
                                  size_tt count,
                                  size_tt placement_mode);
    Point<> initial_position_mode_weighted(const size_tt *closeVert,
                                           const double *closeVertDist,
                                           size_tt count,
                                           size_tt placement_mode);
    Point<> initial_position_barycenter(const size_tt *closeVert,
                                        size_tt count);
    Point<> initial_position_circle(const size_tt *closeVert,
                                    const size_tt *closeVertDist,
                                    size_tt count);
    Point<> initial_position_circle_weighted(const size_tt *closeVert,
                                             const double *closeVertDist,
                                             size_tt count);
    Point<> initial_position_least_squares(const size_tt *closeVert,
                                           const size_tt *closeVertDist,
                                           size_tt count);
    Point<> initial_position_least_squares_weighted(const size_tt *closeVert,
                                                    const double *closeVertDist,
                                                    size_tt count);
    void select_insertion_anchor_subset(std::vector<size_tt> &anchors,
                                        std::vector<size_tt> &anchorDist,
                                        size_tt targetCount);
    void select_insertion_anchor_subset_weighted(std::vector<size_tt> &anchors,
                                                 std::vector<double> &anchorDist,
                                                 size_tt targetCount);
    void compute_weighted_shortest_paths(size_tt root,
                                         std::vector<double> &dist,
                                         double cutoff = std::numeric_limits<double>::infinity()) const;
    void traverse_weighted_shortest_paths(size_tt root,
                                          double cutoff,
                                          size_tt maxSettled,
                                          const std::function<bool(size_tt, double)> &visitor);
    std::vector<size_tt> ordered_vertices_by_metric(const std::vector<double> &dist,
                                                    size_tt root) const;
    size_tt lgkk_round_budget_for_layer(size_tt mishLayer) const;
    bool should_run_multiscale_lgkk(size_tt activeCount,
                                    size_tt mishLayer) const;
    void clear_lgkk_level_cache();
    void build_lgkk_level_cache(size_tt activeCount,
                                size_tt mishLayer);
    void compute_active_shortest_paths(size_tt sourceIndex,
                                       size_tt activeCount,
                                       std::vector<double> &dist,
                                       std::vector<int> *parent = nullptr);
    std::vector<size_tt> lgkk_choose_local_neighbors(size_tt sourceIndex,
                                                     size_tt activeCount) const;
    std::vector<size_tt> lgkk_choose_landmarks(size_tt sourceIndex,
                                               size_tt activeCount) const;
    void lgkk_refine_level(size_tt activeCount,
                           size_tt mishLayer,
                           size_tt &traceRoundInLevel);
    void capture_trace_snapshot(const char *phase,
                                size_tt activeCount,
                                size_tt roundInLevel);
    void trace_begin_level(size_tt activeCount);
    void trace_after_round(size_tt activeCount, size_tt roundInLevel);
    void trace_finalize(size_tt activeCount, size_tt roundInLevel);
    
    // log to the base 2
    int ilog(int n){
        if( n <= 0 ){
            debug("ERROR: ilog is def only for positive numbers");
            throw std::invalid_argument("ilog is defined only for positive integers");
        }
        int k = 0;
        while(n){
            n /= 2;
            k++;
        }
        return k-1;
    }
};

#endif
