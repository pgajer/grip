// Shared support routines for the GRIP layout engine

#include "DrawGraph.h"


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
