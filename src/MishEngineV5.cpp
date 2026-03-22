// mish_engine_v5() and the supporting routines

#include "DrawGraph.h"

/**
 * @brief Implements the Maximal Independent Set Hierarchy (MISH) engine for graph layout.
 *
 * @details This function, mish_engine_v5(), is the core of the graph layout algorithm.
 * It uses a Maximal Independent Set Filtration (MISF) approach to progressively
 * refine the layout of the graph. The algorithm works through different levels
 * of the MISF, applying force-directed layout techniques at each level.
 *
 * Key features:
 * - Supports both arbitrary dimension layouts.
 * - Uses static variables to maintain state across function calls.
 * - Implements a multi-level approach, working from coarse to fine graph representations.
 * - Applies force-directed algorithms (FR_spring_v2 or KK_spring_v4) for layout refinement.
 * - Manages temperature cooling to control vertex movement.
 *
 * @pre The graph structure must be initialized.
 * @pre MISF must be computed (misfSize and misfLevel must be set).
 * @pre Class members (pos, heat, etc.) must be properly initialized.
 *
 * @post Graph layout is updated and refined.
 * @post createList flag is set based on the algorithm's progress.
 *
 * @note This function is designed to be called repeatedly until the desired layout quality is achieved.
 *
 * @param [in,out] displayPar Controls whether the function performs a single iteration (1) or continues until convergence (0).
 *
 * Key member variables used/modified:
 * @li mish: Array representing the MISF structure.
 * @li misfSize: Array storing the size of each MISF level.
 * @li misfLevel: The current level of MISF being processed.
 * @li pos: Array storing vertex positions in 2D or 3D.
 * @li disp: Array storing vertex displacements.
 * @li heat: Array storing temperature values for each vertex.
 * @li nbrs: Array of neighbor lists for each vertex at each MISF level.
 *
 * @see FR_spring_v2() for the Fruchterman-Reingold force-directed algorithm.
 * @see KK_spring_v4() for the Kamada-Kawai force-directed algorithm.
 * @see bfs_me_v4() for the breadth-first search used in initialization.
 * @see update_Local_Temp_v2() for temperature update mechanism.
 *
 * @warning This function modifies global state and should be used carefully in multi-threaded environments.
 *
 * Time Complexity: O(rounds * csize * (V + E)), where V is the number of vertices, E is the number of edges,
 * csize is the current active set size, and rounds is the number of iterations at each level.
 *
 * Space Complexity: O(V), primarily due to the various arrays used to store vertex information.
 *
 * @todo Implemente parallelization for force calculations to improve performance on large graphs.
 * @todo Investigate adaptive schemes for determining optimal number of rounds at each level.
 */
void DrawGraph::mish_engine_v5()
{
    #define DEBUG_mish_engine_v5 0

#if DEBUG_mish_engine_v5
    std::cerr << "Entering mish_engine_v5()  misfLevel: " << misfLevel
              << " displayPar: " << displayPar
              << std::endl;
#endif

    bool firstRound = true;
    size_tt current_size = numOfInitVert; // Current size of the active vertex set
    size_tt round_counter = 0;            // Counter for the current round
    size_tt trace_round_in_level = 0;     // Completed rounds in the current level
    bool loop = true;

#if DEBUG_mish_engine_v5
    int loop_itr = 0; // only for debugging !!!
#endif

    while( loop && prevSize != current_size ){

        // If displayPar is set, only perform one iteration
        if(displayPar)
            loop = false;
#if DEBUG_mish_engine_v5
        std::cerr << "loop_itr: " << loop_itr
                  << " round_counter: " << round_counter
                  << " prevSize: " << prevSize
                  << " current_size: " << current_size
                  << " loop: " << loop
                  << " displayPar: " << displayPar << std::endl;
        loop_itr++;
#endif

        if( firstRound ){
#if DEBUG_mish_engine_v5
            std::cerr << "In firstRound" << std::endl;
#endif
            // Initialize on the first round
            firstRound = false;
            rounds = initRounds;
            prevMishSize = 0;

            // Assign random positions to the initial set of vertices
            for(size_tt i = 0; i < current_size; i++){
                pos[mish[i]] = rand_Point();
                baricenter += pos[mish[i]];
            }
            // Center the layout around the origin
            baricenter /= (coord_t)numOfInitVert;
            for(size_tt i = 0; i < current_size; i++)
                pos[mish[i]] -= baricenter;
            trace_round_in_level = 0;
            trace_begin_level(current_size);

        } else if ( round_counter == rounds ){
            // Prepare for the next MISF level
            round_counter = 0;
            prevSize = current_size;
            prevMishSize = misfSize[misfLevel];

            if( prevMishSize < numOfInitVert )
                prevMishSize = numOfInitVert;

            // Determine the new size for the next level
            if(misfSize[misfLevel] != numOfVert)
                current_size = misfSize[--misfLevel];
            else if(prevSize != numOfVert)
                current_size = numOfVert;
            else
                break;

            // Determine the number of rounds for the new level
            rounds = sched3(current_size, 0, initRounds, numOfVert, finalRounds);

            // Reset initialization constants for vertices
            for(size_tt i = 0; i < prevSize; i++)
                  heat[mish[i]] = tinit;

            // Perform BFS for new vertices
            for(size_tt i = prevSize; i < current_size; i++)
                bfs_me_v4(mish[i]);
            trace_round_in_level = 0;
            trace_begin_level(current_size);
        }

        // Perform force-directed layout if not creating a list and within round limit
        if( !createList && round_counter++ < rounds ){
            for(size_tt i = 0; i < current_size; i++){
                if(!misfLevel)
                    FR_spring_v2(mish[i], nbrs[mish[i]][misfLevel], misfLevel);
                else
                    KK_spring_v4(mish[i], nbrs[mish[i]][misfLevel], misfLevel);

                update_Local_Temp_v2(mish[i]);
                oldDisp[mish[i]] = disp[mish[i]];
                oldDispNorm[mish[i]] = dispNorm[mish[i]];
                disp[mish[i]] *= (coord_t)heat[mish[i]];

                if(dispNorm[mish[i]]){
                    disp[mish[i]] /= dispNorm[mish[i]];
                }
            }
            for(size_tt i = 0; i < current_size; i++)
                pos[mish[i]] += disp[mish[i]];
            trace_round_in_level = round_counter;
            trace_after_round(current_size, trace_round_in_level);
        }
    }

    trace_finalize(current_size, trace_round_in_level);

    // Handle list creation flag
    if( loop && listSwitch ) {
        createList = true;
        listSwitch = false;
    } else
        createList = false;
#if DEBUG_mish_engine_v5
    std::cerr << "Leaving mish_engine_v5()\n";
#endif
}


//****************************************************************
//
//       bfs_me_init_v2
//
//       version of bfs() utilizing nbrs 3D array
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
//       version of bfs() utilizing nbrs 3D array
//
//****************************************************************
void DrawGraph::bfs_me_v4(size_tt root)
{
    size_tt bottomNbrsLayer = 0; // the index i of the first free nbrs[root][i]
    size_tt vert;
    size_tt currentDepth = 0;
    FastQueue<size_tt> vertDepthQueue(4*numOfVert);// an array based queue
    
    size_tt numOfCloseVert = 3;
    std::vector<size_tt> closeVert(numOfCloseVert);
    std::vector<size_tt> closeVertDist(numOfCloseVert);
    size_tt closeVertItr = 0;
    bool closeVertDone = false;

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

                if( !closeVertDone && closeVertItr < numOfCloseVert &&
                    vertDepth[overt] > vertDepth[root]){
                    closeVertDist[closeVertItr] = currentDepth;
                    closeVert[closeVertItr++] = overt;

                    if( closeVertItr == numOfCloseVert ){
                        closeVertDone = true;

                        pos[root] = initial_position(closeVert.data(),
                                                     closeVertDist.data(),
                                                     numOfCloseVert);

                        // set oldDisp vector to the average of
                        // disp vectors of the closest centers
                        oldDisp[root] += oldDisp[closeVert[0]];
                        oldDisp[root] += oldDisp[closeVert[1]];
                        oldDisp[root] += oldDisp[closeVert[2]];
                        oldDisp[root] /= 3;
                        oldDispNorm[root] = oldDisp[root].norm();
                        
                        size_tt itr = 0;
                        while(itr++ < 3){
                            KK_spring_local(root, closeVert.data(),
                                            closeVertDist.data(), numOfCloseVert);
                            update_Local_Temp_v3(root, r, s);
                            oldDisp[root] = disp[root];
                            oldDispNorm[root] = dispNorm[root];
                            disp[root] *= (coord_t)heat[root];
                            if(dispNorm[root])
                                disp[root] /= dispNorm[root];
                            pos[root] += disp[root];
                        }
                    }
                }//end of if( !closeVertDone ...
            }            
        }
    } while ( (!closeVertDone || bottomNbrsLayer <= vertDepth[root]) &&
              !vertDepthQueue.is_empty( ) );
    delete [] nbrCounter;
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
    double norm2;// square of the Eucleadian distance between vert and overt
    
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
//    local beautification used for determining initial positions
//    of "new vertices"
//
//    This version destoys bfs if destroyBFS flag is true.
//
//    the force is calculated using the whole bfs tree, so it is
//    assumed that it is guaged correctly.
//    bfs (=bfsVectQueue) is shifted by 'shift', so one needs to
//    take this into account when computing the force vector
//
//**************************************************************
void DrawGraph::KK_spring_local(const size_tt vert,
                                size_tt *closeVert,
                                size_tt *closeVertDist,
                                size_tt size)
{
    coord_t norm2;   // its norm squared
//    queue<size_tt> currentQueue;

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
        coord_t r = 0.15;
        coord_t s = 3.0;
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
//        r = 0.15;
//        s = 3.0;        
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
//    double dist2;// square of the graph theoretic dist between vert and overt
    double norm2;// square of the Eucleadian distance between vert and overt
    
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
