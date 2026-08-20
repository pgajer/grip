// MishEngineV5.cpp - mish_engine_v5() and the supporting routines

#include "DrawGraph.hpp"

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
    size_tt *nbrCounter = new size_tt[vertDepth[root]];
    nbrs[root] = new size_tt*[vertDepth[root]+2];
    for(size_tt i=0; i <= vertDepth[root]; i++){
        // 1 slot for vertex and 1 (the next one)
        // for its distance from the root
        nbrs[root][i] = new size_tt[2*nbr[i]];
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
                    i <= min(vertDepth[overt], vertDepth[root]); i++){
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
    size_tt closeVert[numOfCloseVert];
    size_tt closeVertDist[numOfCloseVert];
    size_tt closeVertItr = 0;
    bool closeVertDone = false;

    marked[root] = root;
    vertDepthQueue.enqueue(root);
    vertDepthQueue.enqueue(currentDepth);
    
    // memory allocation for rootNbrs
    nbrs[root] = new size_tt*[vertDepth[root]+1];
    size_tt *nbrCounter = new size_tt[vertDepth[root]+1];
    for(size_tt i=0; i <= vertDepth[root]; i++){
        // 1 slot for vertex and 1 (the next one)
        // for its distance from the root
        nbrs[root][i] = new size_tt[2*nbr[i]];
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
                    i <= min(vertDepth[overt], vertDepth[root]); i++)
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

                        if( dim == 4 ){
                            pos4D[root] = pos4D[closeVert[0]];
                            pos4D[root] += pos4D[closeVert[1]];
                            pos4D[root] += pos4D[closeVert[2]];
                            pos4D[root] /= 3;

                            // set oldDisp vector to the average of
                            // disp vectors of the closest centers
                            oldDisp4D[root] += oldDisp4D[closeVert[0]];
                            oldDisp4D[root] += oldDisp4D[closeVert[1]];
                            oldDisp4D[root] += oldDisp4D[closeVert[2]];
                            oldDisp4D[root] /= 3;
                            oldDispNorm[root] = oldDisp4D[root].norm();
                            
                            size_tt itr = 0;
                            while(itr++ < 3){
                                KK_spring_local(root, closeVert,
                                                closeVertDist, numOfCloseVert);
                                update_Local_Temp_v3(root, r, s);
                                oldDisp4D[root] = disp4D[root];
                                oldDispNorm[root] = dispNorm[root];
                                disp4D[root] *= (coord_t)heat[root];
                                if(dispNorm[root])
                                    disp4D[root] /= dispNorm[root];
                                pos4D[root] += disp4D[root];
                            }
                        } else {
                            pos[root] = pos[closeVert[0]];
                            pos[root] += pos[closeVert[1]];
                            pos[root] += pos[closeVert[2]];
                            pos[root] /= 3;

                            // set oldDisp vector to the average of
                            // disp vectors of the closest centers
                            oldDisp[root] += oldDisp[closeVert[0]];
                            oldDisp[root] += oldDisp[closeVert[1]];
                            oldDisp[root] += oldDisp[closeVert[2]];
                            oldDisp[root] /= 3;
                            oldDispNorm[root] = oldDisp[root].norm();
                            
                            size_tt itr = 0;
                            while(itr++ < 3){
                                KK_spring_local(root, closeVert,
                                                closeVertDist, numOfCloseVert);
                                update_Local_Temp_v3(root, r, s);
                                oldDisp[root] = disp[root];
                                oldDispNorm[root] = dispNorm[root];
                                disp[root] *= (coord_t)heat[root];
                                if(dispNorm[root])
                                    disp[root] /= dispNorm[root];
                                pos[root] += disp[root];
                            }
                        }
                    }
                }//end of if( !closeVertDone ...
            }            
        }
    } while ( (!closeVertDone || bottomNbrsLayer <= vertDepth[root]) &&
              !vertDepthQueue.is_empty( ) );
}


//****************************************************************
//
//     bfs_cmisf()
// 
//     version of bfs_cm4 (used by create_misf()), which
//     rewritten to use misfSize array
//
//****************************************************************
void DrawGraph::bfs_cmisf(size_tt root,
                          queue<size_tt> *bfsVectQueue,
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

//**************************************************************
//
//	mish_engine_v5()
//
//      "maximal independent set filtration (former hierarchy)" engine
//
//      main changes:
//      utilizes bfs_me_v3(), KK_spring_v4(), nbrsArr[], vertDepth[]
//
//**************************************************************
void DrawGraph::mish_engine_v5()
{
    static bool firstRound = true;
    static size_tt csize = numOfInitVert;
    static size_tt ctr = 0;
    bool loop = true;
    
    while( loop && prevSize != csize ){
        if(displayPar)
            loop = false;

        if( firstRound ){
            firstRound = false;
            rounds = initRounds; 
            prevMishSize = 0;
            
            // Assign random positions to the initial set of vertices
            if( dim == 4 ){
                for(size_tt i = 0; i < csize; i++){
                    pos4D[mish[i]] = rand_Point4D();
                    baricenter4D += pos4D[mish[i]];
                }
            } else {
                for(size_tt i = 0; i < csize; i++){
                    pos[mish[i]] = rand_Point();
                    baricenter += pos[mish[i]];
                }
            }
            if( dim == 4 ){
                baricenter4D /= (coord_t)numOfInitVert;
                for(size_tt i = 0; i < csize; i++)
                    pos4D[mish[i]] -= baricenter4D;
            } else {
                baricenter /= (coord_t)numOfInitVert;
                for(size_tt i = 0; i < csize; i++)
                    pos[mish[i]] -= baricenter;
            }
            
            size_tt numOfNbrsLim = numOfInitVert;
            debug("csize="<< numOfInitVert
                  << ", nbrsLim="<< numOfNbrsLim
                  << ", rounds=" << rounds);

        } else if ( ctr == rounds ){
            ctr = 0;
            prevSize = csize;
            prevMishSize = misfSize[misfLevel];

            if( prevMishSize < numOfInitVert )
                prevMishSize = numOfInitVert;
            if(misfSize[misfLevel] != numOfVert)
                csize = misfSize[--misfLevel];
            else if(prevSize != numOfVert)
                csize = numOfVert;
            else
                break;

            // DETERMINING rounds
            rounds = sched3(csize,
                            0, initRounds,
                            numOfVert, finalRounds);        


            // RESETTING INITIALIZATION CONSTANTS FOR NEW AND OLD
            // VERTICES
            for(size_tt i = 0; i < prevSize; i++)
                  heat[mish[i]] = tinit;

            for(size_tt i = prevSize; i < csize; i++)
                bfs_me_v4(mish[i]);
            
            debug("csize="<< csize
                  <<", nbr["<<misfLevel<<"]="<< nbr[misfLevel]
                  <<", rounds=" << rounds);

        }// end of PREPROCESSING SECTION if( firstRound ){} else if( ...

    
        if( !createList && ctr++ < rounds ){
            if( dim == 4 ){
                for(size_tt i = 0; i < csize; i++){
                    //perform local force-directed modifications
                    // rounds-ctr is used here as a destroyFlag
                    // not used in this version of KK_spring
                    if(!misfLevel)
                        FR_spring_v2(mish[i], nbrs[mish[i]][misfLevel], misfLevel);
                    else
                        KK_spring_v4(mish[i], nbrs[mish[i]][misfLevel], misfLevel);
                    update_Local_Temp_v2(mish[i]);
                    oldDisp4D[mish[i]] = disp4D[mish[i]];
                    oldDispNorm[mish[i]] = dispNorm[mish[i]];
                    disp4D[mish[i]] *= (coord_t)heat[mish[i]];
                    
                    if(dispNorm[mish[i]]){
                        disp4D[mish[i]] /= dispNorm[mish[i]];
                    }
                }
                for(size_tt i = 0; i < csize; i++)
                    pos4D[mish[i]] += disp4D[mish[i]];

            } else {
                for(size_tt i = 0; i < csize; i++){
                    //perform local force-directed modifications
                    // rounds-ctr is used here as a destroyFlag
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
                for(size_tt i = 0; i < csize; i++)
                    pos[mish[i]] += disp[mish[i]];
            }// end of if( dim == 4 )
        }// end of if( !createList
    } // end of while( loop && prevSize != ...
    
    if( loop && listSwitch ) {
        createList = true;
        listSwitch = false;
    } else
        createList = false;

//    debug("leaving mish5");
}

//**************************************************************
//
//     create_misf()
//
//     building a maximal independent set filtration
//     is it version 5 of create_mish()
//     it creates also an array vertDepth
//
//**************************************************************
void DrawGraph::create_misf(){

    size_tt depthLim;
    size_tt shift=0;
    queue<size_tt> **bfsVectQueue = new queue<size_tt>*[numOfVert];

    misfLevel = 1;
    misfSize[misfLevel-1] = numOfVert;
    size_tt mishSizeCurrLevel = numOfVert;
 
    size_tt itr;
    do{
        size_tt mishSizePrevLevel = mishSizeCurrLevel;
        mishSizeCurrLevel = 0;
        itr = 0;     // independent set iterator (runs through) mish

        // the following while loop creates elements of mish[misfLevel]
        while( mishSizePrevLevel > numOfInitVert &&
               itr < mishSizePrevLevel ){
            size_tt vert = itr + graph.fast_Rand() % (mishSizePrevLevel - itr);
            
            swap(mish[vert], mish[itr]);
            inv[mish[vert]] = vert;
            inv[mish[itr]] = itr;
            
            swap(mish[mishSizeCurrLevel], mish[itr]);
            inv[mish[mishSizeCurrLevel]] = mishSizeCurrLevel;
            inv[mish[itr]] = itr;
            itr++;

            // this is a new element of mish[misfLevel]
            size_tt newEl = mish[mishSizeCurrLevel++];

            vertDepth[newEl] = misfLevel;
            
            // computing maximal independent set
            if( misfLevel == 1 ){
                size_tt adj;
                for(size_tt i=0; i<deg[newEl] && itr<numOfVert; i++){
                    adj = graph.adjList[newEl+1][i];
                    if(inv[adj] >= itr){
                        swap(mish[itr], mish[inv[adj]]);
                        inv[mish[inv[adj]]] = inv[adj];
                        inv[mish[itr]] = itr;
                        itr++;
                    }
                }
            } else {
                for(size_tt shiftedDepth=0;
                    shiftedDepth < depthLim-shift+1;
                    shiftedDepth++){
                    while( !bfsVectQueue[newEl][shiftedDepth].empty() ){
                        size_tt adj = bfsVectQueue[newEl][shiftedDepth].front();
                        bfsVectQueue[newEl][shiftedDepth].pop();
                        if(inv[adj] >= itr){
                            swap(mish[itr], mish[inv[adj]]);
                            inv[mish[inv[adj]]] = inv[adj];
                            inv[mish[itr]] = itr;
                            itr++;
                        }
                    } 
                }// end of for(size_tt shiftedDepth=0 ...
            }// end of if( mishSize == 1 ) .. else ...
        }// end of while( mishSize[misfLevel-1] > numOfInitVert && ...


        depthLim = (size_tt)pow(2, misfLevel);
        depthLim = min(depthLim, (size_tt)(numOfVert-1));

        assert(misfLevel < log_2_n);
        misfSize[misfLevel] = mishSizeCurrLevel;

        // creating bfsVectQueue structores for the elements
        // of mish[misfLevel]

        for(int i = 0; i < numOfVert; i++) 
            marked[i] = -1;
        shift = (size_tt)pow(2, misfLevel-1)+1;
        for(size_tt i=0; i < mishSizeCurrLevel; i++){
            bfsVectQueue[mish[i]] = new queue<size_tt>[depthLim-shift+1];
            bfs_cmisf(mish[i], bfsVectQueue[mish[i]], shift, depthLim);
        }
        misfLevel++;
    } while(itr);

    assert(misfLevel <= log_2_n);
    // set mishLevel to the first level of mish which is >= numOfInitVert
    misfLevel -= 1;
    while( misfSize[misfLevel] < numOfInitVert )
        misfLevel--;

    size_tt v=0;
    if(misfSize[misfLevel] > numOfInitVert){
        while( v < numOfInitVert )
            vertDepth[mish[v++]] = misfLevel+1;
        misfLevel++;
        misfSize[misfLevel] = numOfInitVert;
    }
    
    debug("misfLevel=" << misfLevel);
    initMishHeight = misfLevel;

    for(int i=0; i<numOfVert;i++)
        delete [] bfsVectQueue[i];
    delete [] bfsVectQueue;
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

    if( dim == 4 ){
        disp4D[vert].set_to_zero();
        for(size_tt i = 0; i < 2*nbr[misfLayer]; i += 2){
            overt = *ptr++;
            dist2 = (double)(*ptr) * (*ptr);
            ptr++;
            if(!dist2)
                continue;
            vect4D.set_to_zero();
            vect4D += pos4D[overt];
            vect4D -= pos4D[vert];
            norm2 = (double)vect4D.fnorm2();
            vect4D *= (float)(norm2/(dist2 * edge2) - 1);
            disp4D[vert] += vect4D;
        }
        
        float norm = disp4D[vert].fnorm();
        dispNorm[vert] = ROUND_L(norm);
        
        if(dispNorm[vert]){
            disp4D[vert] *= edge/norm;    
            dispNorm[vert] = disp4D[vert].norm();
        }
    } else {
        disp[vert].set_to_zero();
        for(size_tt i = 0; i < 2*nbr[misfLayer]; i += 2){
            overt = *ptr++;
            dist2 = (double)(*ptr) * (*ptr);
            ptr++;
            if(!dist2)
                continue;
            vect.set_to_zero();
            vect += pos[overt];
            vect -= pos[vert];
            norm2 = (double)vect.fnorm2();
            vect *= (float)(norm2/(dist2 * edge2) - 1);
            disp[vert] += vect;
        }
        
        float norm = disp[vert].fnorm();
        dispNorm[vert] = ROUND_L(norm);
        
        if(dispNorm[vert]){
            disp[vert] *= edge/norm;    
            dispNorm[vert] = disp[vert].norm();
        }
    }
}

//**************************************************************
//
//    KK_spring_local()
//
//    local beautification used for determining initial positions
//    of "new vertices"

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

    if( dim == 4 ){
        disp4D[vert].set_to_zero();
        for(size_tt i = 0; i < size; i++){
            vect4D.set_to_zero();
            vect4D = pos4D[closeVert[i]] - pos4D[vert];
            norm2 = vect4D.norm2();
            vect4D *=
                ((float)norm2/(closeVertDist[i]*closeVertDist[i]*edge2)-1);
            disp4D[vert] += vect4D;
        }
        dispNorm[vert] = disp4D[vert].norm();
    } else {
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
    size_tt temp = heat[vert];
    unsigned long normOldDisp = (unsigned long)oldDispNorm[vert];
    unsigned long normNewDisp = (unsigned long)dispNorm[vert]; 
    
    if( normOldDisp != 0 && normNewDisp != 0 ){
        coord_t scalProd;
        if( dim == 4 )
            scalProd = disp4D[vert] * oldDisp4D[vert];
        else
            scalProd = disp[vert] * oldDisp[vert];
        float cos = (float)scalProd/(normOldDisp * normNewDisp);
        float r = 0.15;
        float s = 3;
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
void DrawGraph::update_Local_Temp_v3( size_tt vert, float r, float s)
{
    size_tt temp = heat[vert];
    unsigned long normOldDisp = (unsigned long)oldDispNorm[vert];
    unsigned long normNewDisp = (unsigned long)dispNorm[vert]; 
    
    if( normOldDisp != 0 && normNewDisp != 0 ){
        coord_t scalProd;
        if( dim == 4 )
            scalProd = disp4D[vert] * oldDisp4D[vert];
        else
            scalProd = disp[vert] * oldDisp[vert];
        float cos = (float)scalProd/(normOldDisp * normNewDisp);
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

    if( dim == 4 ){
        disp4D[vert].set_to_zero();

        // attractive force calculation
        size_tt deg = graph.adjList[0][vert];
        for (size_tt adjVert = 0; adjVert < deg; adjVert++){
            overt = graph.adjList[vert+1][adjVert];
            vect4D.set_to_zero();
            vect4D += pos4D[overt];
            vect4D -= pos4D[vert];
            norm2 = (double)vect4D.fnorm2();
            vect4D *= (float)vect4D.norm2();
            vect4D /= edge2;
            disp4D[vert] += vect4D;
        }

        // repulsive force calculation
        size_tt locNbr = 2*nbr[misfLayer];
        for(size_tt i = 0; i < locNbr; i += 2){
            overt = *ptr++;
            vect4D.set_to_zero();
            vect4D += pos4D[vert];
            vect4D -= pos4D[overt];
            norm2 = (double)vect4D.fnorm2();
            if(!norm2)
                continue;
            vect4D *= (float)(fedge2/norm2);
            disp4D[vert] += vect4D;
        }
        float norm = disp4D[vert].fnorm();
        dispNorm[vert] = ROUND_L(norm);

        if(dispNorm[vert]){
            disp4D[vert] *= edge/norm;    
            dispNorm[vert] = disp4D[vert].norm();
        }
    } else {
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
            vect /= edge2;
            disp[vert] += vect;
        }

        // repulsive force calculation
        size_tt locNbr = 2*nbr[misfLayer];
        for(size_tt i = 0; i < locNbr; i += 2){
            overt = *ptr++;
            vect.set_to_zero();
            vect += pos[vert];
            vect -= pos[overt];
            norm2 = (double)vect.fnorm2();
            if(!norm2)
                continue;
            vect *= (float)(fedge2/norm2);
            disp[vert] += vect;
        }
        float norm = disp[vert].fnorm();
        dispNorm[vert] = ROUND_L(norm);

        if(dispNorm[vert]){
            disp[vert] *= edge/norm;    
            dispNorm[vert] = disp[vert].norm();
        }
    }
}
