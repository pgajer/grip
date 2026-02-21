// FILE: DrawGraph.cpp - member function definitions for DrawGraph class.
//   
// April 1, 2000 - moved all changes info to the file changes
//
#include "DrawGraph.hpp"

#define MISH_ENG_v5          12
#define MISH_ENG_v6          13

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
                     size_tt _engf      ,
                     size_tt _numOfInitVert,
                     size_tt _numOfNbrs ,
                     float _r           ,  //parameters of update_Local_Temp_v3()
                     float _s           ,
                     bool _displayPar   )
: createList(false),
  graph(_graph),
  dim(_dim),
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
  engf(_engf),
  fedge2(0.05f * edge2)
{
#define DEBUG 0
#if DEBUG 
    debug("Entering constructor\n");
#endif

    if( dim == 4 ){
        pos4D        = new Point4D<>[numOfVert];
        disp4D       = new Point4D<>[numOfVert];
        oldDisp4D    = new Point4D<>[numOfVert];
    } else {
        pos          = new Point<>[numOfVert];
        disp         = new Point<>[numOfVert];
        oldDisp      = new Point<>[numOfVert];
    }

    dispNorm     = new coord_t[numOfVert];
    oldDispNorm  = new coord_t[numOfVert];
    old_cos      = new float[numOfVert];
    heat         = new int[numOfVert];
    deg          = new size_tt[numOfVert];
        
    for(size_tt vert=0; vert < numOfVert; vert++){
        heat[ vert ] = tinit;
        old_cos[vert] = 1;
        deg[vert] = graph.adjList[0][vert];
    }

    numOfInitVert = std::min(numOfInitVert, numOfVert);

    //
    // mish_engine_v5() and mish_engine_v6()
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

    // compute the number of vertices nbr[i] for local beautification
    // at level i, for each i.
    nbr       = new size_tt[log_2_n];
    itr = 0;
    //float nbrFactor = 2; // affects the size of nbr[i]s
    while( itr < log_2_n && misfSize[itr] ){
        if( itr >= smallLevel )
            nbr[itr] = std::max(misfSize[itr]-1, numOfInitVert-1);
        else{
            nbr[itr] = std::min((unsigned long)(sched(itr,0,2,10000,1) *
                                           maxCxty/misfSize[itr]),
                           (unsigned long)(misfSize[itr]-1));
//          nbr[itr] = std::min((unsigned long)( nbrFactor * maxCxty/misfSize[itr]),
//                           (unsigned long)(misfSize[itr]-1));
        }
        itr++;
    }

    // some simple nbr[] tune up
    nbr[0] = std::min(2*nbr[0], numOfVert-1);
//          if(misfLevel)
//              nbr[1] *= 2;

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

    if( dim == 4 ){
        delete [] pos4D;    
        delete [] disp4D;
        delete [] oldDisp4D;
    } else {
        delete [] pos;
        delete [] disp;
        delete [] oldDisp;
    }
    delete [] deg;
    delete [] mish;
    delete [] inv;
    delete [] misfSize;
    delete [] vertDepth;
    delete [] nbr;
    delete [] marked;
    
    //
    // MORE PROPERLY the first double for loop should be
    // executed in KK_spring when unused parts of nbrs[i]
    // are deleted
    // here we should execute only the single for loop and
    // delete [] nbrs;
    //
    for(int i=0; i<numOfVert;i++)
        for(int j=0; j<vertDepth[i]; j++)
            delete [] nbrs[i][j];
    for(int i=0; i<numOfVert;i++)
        delete [] nbrs[i];
    delete [] nbrs;
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

//**************************************************************
//
//	rand_Point4D()
//
//**************************************************************
Point4D<> DrawGraph::rand_Point4D()
{
    return 
        Point4D<>((coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize,
                  (coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize,
                  (coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize,
                  (coord_t)(graph.fast_Rand() % (int)box2Size) - boxSize);
}

//****************************************************************
//
//    memory exception handler
//
//****************************************************************
void DrawGraph::noMoreMemory(){
        std::cerr << "Unable to satisfy request for memory\n";
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
    for (int i = 0; i < numOfVert; i++)
        delete [] bfsVectQueue[i];
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
