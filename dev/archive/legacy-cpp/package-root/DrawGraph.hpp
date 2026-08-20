// DrawGraph.hpp a header file for DrawGraph class

#ifndef DRAW_GRAPH_HPP
#define DRAW_GRAPH_HPP

#include <iostream>
#include <new>
#include <cmath>
#include <queue>

#include "Point.hpp"
#include "Point4D.hpp"
#include "Graph.hpp"
#include "FastQueue.hpp"
#include "Debug.hpp"

typedef unsigned short size_tt;

#define ROUND(a)       ((a)>0 ? (int)((a)+0.5) : -(int)(0.5-(a)))
#define ROUND_L(a) ((a)>0 ? (unsigned long)((a)+0.5) : -(unsigned long)(0.5-(a)))

//**************************************************************
//
//	Class name : DrawGraph - compact version
//      only with mish_engine versions 5 and 6
//
//**************************************************************
class DrawGraph 
{
  public:
      friend class MesaPlot;
    
    DrawGraph(const Graph &_graph,
              size_tt _dim         ,
              size_tt _rounds      ,
              size_tt _finalRounds ,
              size_tt _tinit_factor,
              size_tt _engf        ,
              size_tt _numOfInitVert,
              size_tt _numOfNbrs   ,
              float _r             ,  //parameters of update_Local_Temp_v3()
              float _s             ,
              bool _displayPar     );
    
    ~DrawGraph();

    //
    // ENGINES
    //
    //----- mish_engine_v5 --------------//
    void mish_engine_v5();
    size_tt bfs_me_init_v2(size_tt root);
    void bfs_me_v4(size_tt root);
    void bfs_cmisf(size_tt root,
                   std::queue<size_tt> *bfsVectQueue,
                   size_tt shift,   // bfsVector starts at dist shift from root
                   size_tt depthLim);
    void create_misf();  // version of create_mish() creating also vertDepth
    void KK_spring_v4(const size_tt root,
                      size_tt *rootNbrsLayer,
                      size_tt mishLayer);
    void FR_spring_v2(const size_tt vert,
                      size_tt *vertNbrs,
                      size_tt misfLayer);
    void KK_spring_local(const size_tt vert,
                         size_tt *closeVert,
                         size_tt *closeVertDist,
                         size_tt size);
    void update_Local_Temp_v2( size_tt vert );
    void update_Local_Temp_v3( size_tt vert, float r, float s );
    
    //----- mish_engine_v6 --------------//    
    void mish_engine_v6(); // version utilizing repulsive forces
    // FR_spring() is used in mish_engine_v6() and utilizes
    // attractive/repulsive force schedule
    void FR_spring(const size_tt root,
                   size_tt *rootNbrsLayer,
                   size_tt mishLayer);
    

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
    Point4D<> *get_Pos4D() const { return pos4D; }
    coord_t get_Edge(){ return edge; }
    coord_t get_Edge2(){ return edge2; }
    size_tt get_inv(size_tt vert){return inv[vert];}
    size_tt get_prevMishSize(){return prevMishSize;}
    coord_t dist(const Point<> &p, const Point<> &q);
    
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

    const size_tt edge;               // desired edge length [128]
    const size_tt edge2;              // = edge * edge
    size_tt tinit;              // desired minimal temp [3]
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
    
    Point4D<> center4D;      // = \sum pos[vert]
    Point4D<> baricenter4D;  // center / numOfVert
    Point4D<> *pos4D;        // positions of vertices vector
    Point4D<> *disp4D;       // displacement vector
    Point4D<> *oldDisp4D; 

    size_tt *deg;             // array of degrees
    float *old_cos;
    int *maxNorm;              // vector of maximal norms of disp
    int *heat;                   // local temperature

    //
    // variable used specifically in mish_engines
    //
    size_tt *mish;               // maximal independent set hierarchy
    size_tt *inv;                // "inverse" of mish
    size_tt *misfSize;  // mish_engine_v5 analogue of mishSize
    size_tt log_2_n;    // log_2(n)
    size_tt *vertDepth; // array of depths of vertices
    float AvgDeg;
    unsigned long maxCxty;
    unsigned long initCxty;
    size_tt smallLevel;
    size_tt *nbr;       // array of num of nbrs for each level of misf
    size_tt ***nbrs;    // storage of nbrs for corresponding levels
    size_tt prevSize;
    
    size_tt prevMishSize;
    size_tt diam;                         // graph's diameter
    size_tt numOfInitVert;
    size_tt numOfNbrs;
    float r, s;                           // parameters of update_Local_Temp_v3()
    bool listSwitch;
    bool displayPar;
    size_tt engf;
    float fedge2;
    size_tt misfLevel;    
    size_tt initMishHeight;
    Point<> vect;    // an auxiliary vector holding vert1 - vert2
    Point4D<> vect4D;    // and its 4D analogue

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

    size_t norm2(const Point<> &pt1, const Point<> &pt2){
        return ((pt1.getX() - pt2.getX()) * (pt1.getX() - pt2.getX()) +
                (pt1.getY() - pt2.getY()) * (pt1.getY() - pt2.getY()) +
                (pt1.getZ() - pt2.getZ()) * (pt1.getZ() - pt2.getZ()));
    }

    size_t norm2(size_tt u, size_tt v){
        return norm2(pos[u], pos[v]);
    }
    
    size_tt get_BoxSize(){ return boxSize; }

    float sign( coord_t x )
        { 
        if( x > 0 )
            return 1.0;
        else if( x< 0 )
            return -1.0;
        else
            return 0.0;
        }

    float sign( float x )
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
    Point4D<> get_Baricenter4D()
        {
            Point4D<> center4D;
            for(size_tt vert = 0; vert < numOfVert; vert++)
                center4D += pos4D[vert];
            return center4D /(coord_t)numOfVert;
        }

    // random point
    Point<> rand_Point();
    Point4D<> rand_Point4D();
    
    // log to the base 2
    int ilog(int n){
        if( n <= 0 ){
            debug("ERROR: ilog is def only for positive numbers");
            exit(1);
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

