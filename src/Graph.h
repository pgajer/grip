// Graph class declaration.

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cstdlib>
#include <ctime>
#include <vector>
#include <utility>
#include <cstdint>

#include "Point.h"

using size_tt = uint32_t;

//**************************************************************
//
//	Class name: Graph
//
//	Description: A class for generating graphs and performing
//      some operations on them.
//
//**************************************************************
class Graph
{
public:
    friend class DrawGraph;
    friend class MesaPlot;

    Graph() : numOfVert(0), adjList(nullptr), adjWeight(nullptr) {}
    ~Graph(){ clear(); }
    void rand_Graph( size_tt h, size_tt w );
    void complete_Graph( size_tt numOfVert );
    void path_Graph( size_tt numOfVert );
    void cycle_Graph( size_tt numOfVert );

    void double_Graph();
    void pow2( size_tt exp ); // doubling the graph exp times
    void hyper_Cube( size_tt dim );
    void square_Cylinder( size_tt h, size_tt w );
    void square_Torus( size_tt height);
    void mesh(size_tt h);
    void meshX(size_tt h);
    void addEdge(size_tt a, size_tt b);
    size_tt sierpinski_recurse3D(int maxLevel, int currentLevel,
                                size_tt a, size_tt b, size_tt c,
                                size_tt d, bool flag);
    void sierpinski_recurse(int maxLevel, int currentLevel,
                           size_tt a, size_tt b, size_tt c);
    void sierpinski(size_tt h, size_tt density);
    void meshT(size_tt h);  
    void meshTX(size_tt h);  
    void mesh_Graph(size_tt h, size_tt w);
    void torus( size_tt h, size_tt w);
    void twistedTorus( size_tt h, size_tt w,
                       size_tt t1, size_tt t2);
    void torus_Graph(size_tt h, size_tt w);
    void tree(size_tt h, size_tt w);
  
    void from_edge_list(size_tt n,
                        const std::vector<std::pair<size_tt, size_tt>> &edges,
                        const std::vector<coord_t> *weights = nullptr);
    void from_adj_list(size_tt n,
                       const std::vector<std::vector<size_tt>> &adj,
                       const std::vector<std::vector<coord_t>> *weights = nullptr);
    void clear();
    
    size_tt **get_adjList() const { return adjList; }
    
    size_tt get_numOfVert() const { return numOfVert; }
    size_tt get_Deg(size_tt vert) const { return adjList[0][vert]; }
    size_tt get_adjVert(size_tt vert, size_tt adjVert) const
        { return adjList[vert+1][adjVert]; }
    coord_t get_edge_weight(size_tt vert, size_tt adjIndex) const {
        if(!adjWeight) return 1.0;
        return adjWeight[vert+1][adjIndex];
    }
    bool has_weights() const { return adjWeight != nullptr; }
    
    // Return a maximal-degree vertex and report the degree range.
    size_tt get_Max_Deg_Vert(size_tt &max_deg, size_tt &min_deg) const;
    
    // fast random number generator from "Numerical Recipes in C"
    // section "An even quicker generator" p. 284.
    unsigned long fast_Rand() const;
    void sfast_Rand(unsigned int seed) const;
    
    size_tt **rand_cpt_Graph( size_tt numOfVert );
    
    size_tt get_numOfEdges(){
        size_tt numOfEdges = 0;
        for(size_tt vert = 0; vert < numOfVert; vert++)
            numOfEdges += adjList[0][vert];
        return numOfEdges /= 2;
    }
    
private:
    size_tt numOfVert;
    size_tt **adjList; // adjacencey list of the graph
    coord_t **adjWeight; // parallel adjacency weights (optional)


    
    // Swap two integers, coord_t values, or Point<> objects.
    void swap(int &a, int &b);
    void swap(size_tt &a, size_tt &b);    
    void swap(Point<> &a, Point<> &b); 

    // perform a permutation of 'array' of 'len' elements
    // (We can set 'len' to be any number <= length(array).)
    // Put the first newLen elements of 'array' into newArray.
    void rand_Perm(size_tt *array,
                   size_tt len,
                   size_tt *newArray,
                   size_tt newLen);
    
    int min(int a, int b) const { return (a >= b) ? b : a; }

};

#endif

    
