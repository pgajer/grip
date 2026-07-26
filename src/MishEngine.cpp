// GRIP layout engine.

#include "DrawGraph.h"

//**************************************************************
//
//	mish_engine()
//
//      Run multiscale refinement with attractive and repulsive forces.
//
//**************************************************************
void DrawGraph::mish_engine()
{
    if(weightedCore){
        mish_engine_weighted();
        return;
    }

    bool firstRound = true;
    size_tt csize = numOfInitVert;
    size_tt ctr = 0;
    size_tt trace_round_in_level = 0;
    bool loop = true;
    
    while( loop && prevSize != csize ){
        if(displayPar)
            loop = false;

        if( firstRound ){
            firstRound = false;
            rounds = initRounds; 
            prevMishSize = 0;
            currentRoundInLevel = 0;
            
            // Assign random positions to the initial set of vertices
            for(size_tt i = 0; i < csize; i++){
                pos[mish[i]] = rand_Point();
                center += pos[mish[i]];
            }
            baricenter = center /(coord_t)numOfInitVert;
            for(size_tt i = 0; i < csize; i++)
                pos[mish[i]] -= baricenter;
            trace_round_in_level = 0;
            if(!misfLevel && csize == numOfVert){
                for(size_tt i = 0; i < csize; i++)
                    finalAnchorPos[mish[i]] = pos[mish[i]];
                finalAnchorReady = true;
            }
            trace_begin_level(csize);
            
        } else if ( ctr == rounds ){
            if(!createList)
                lgkk_refine_level(csize, misfLevel, trace_round_in_level);
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

            // Determine the refinement budget for the next level.
            rounds = sched3(csize,
                            0, initRounds,
                            numOfVert, finalRounds);        


            // Reset the temperature for previously placed vertices.
            for(size_tt i = 0; i < prevSize; i++)
                  heat[mish[i]] = tinit;
        
            for(size_tt i = prevSize; i < csize; i++)
                bfs_me_v4(mish[i]);
            currentRoundInLevel = 0;
            if(!misfLevel && csize == numOfVert){
                for(size_tt i = 0; i < csize; i++)
                    finalAnchorPos[mish[i]] = pos[mish[i]];
                finalAnchorReady = true;
            }
            trace_round_in_level = 0;
            trace_begin_level(csize);
            
        }// end of PREPROCESSING SECTION if( firstRound ){} else if( ...

    
        if( !createList && ctr++ < rounds ){
            activeVertCount = csize;
            currentRoundInLevel = ctr;
            for(size_tt i = 0; i < csize; i++){
                //perform local force-directed modifications
                // rounds-ctr is used here as a destroyFlag
                if(!misfLevel && finalStageMode == FINAL_STAGE_KK_REPULSE)
                    KK_spring_final(mish[i], nbrs[mish[i]][misfLevel], misfLevel);
                else if(!misfLevel)
                    FR_spring(mish[i], nbrs[mish[i]][misfLevel], misfLevel);
                else
                    KK_spring_v4(mish[i], nbrs[mish[i]][misfLevel], misfLevel);
                
                update_Local_Temp_v2(mish[i]);
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
            }
            for(size_tt i = 0; i < csize; i++)
                pos[mish[i]] += disp[mish[i]];
            trace_round_in_level = ctr;
            trace_after_round(csize, trace_round_in_level);
        }// end of if( !createList
    } // end of while( loop && prevSize != ...

    trace_finalize(csize, trace_round_in_level);
    
    if( loop && listSwitch ) {
        createList = true;
        listSwitch = false;
    } else
        createList = false;
}

//**************************************************************
//
//    FR_spring()
//
//    force calculation
//
//**************************************************************
void DrawGraph::FR_spring(const size_tt vert,
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

    if(!misfLayer)
        add_final_anchor_force(vert);
    
    coord_t norm = disp[vert].fnorm();
    dispNorm[vert] = ROUND_L(norm);

    if(dispNorm[vert]){
        disp[vert] *= edge/norm;    
        dispNorm[vert] = disp[vert].norm();
    }
}
