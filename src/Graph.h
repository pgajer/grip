// Graph.hpp a header file for Graph class

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <utility>
#include <cstdint>

#include "Point.h"
#include "Debug.h"

using size_tt = uint32_t;

#define DEBUG_GRAPH 0

//**************************************************************
//
//	Class name: Graph
//
//	Description: A class for generating graphs and performing
//      some operations on them.
//
//      Date: 9/6/99
//
//      10/23/99: replacing compact adjList by fat adjList, where
//      compact adjList is an adjacency list in which the entries
//      for the i-th vertex contain references only to vertices of
//      index > i, and fat adjacency list is the standard one.
//
//      Author: Pawel Gajer
//
//**************************************************************
class Graph
{
    friend std::ostream &operator<<(std::ostream &output,
                                    const Graph &graph);
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
  
    void print_Cpt_Adj_List_To_File();// print compact Adj list ...
    void read_Cpt_Adj_List_From_File(const char *file);
    //    void read_Graph_From_IWB_File_v1(const char *file);
    void read_Graph_From_IWB_File(const char *file);    
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
    
    // return maximal degree vertex, its degree, and also
    // minimal degree of the graph
    size_tt get_Max_Deg_Vert(size_tt &max_deg, size_tt &min_deg) const;
    
    // fast random number generator from "Numerical Recipes in C"
    // section "An even quicker generator" p. 284.
    //    unsigned long idum;
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


    
    // swap two intergers, coord_t or two Point<>
    void swap(int &a, int &b);
    void swap(size_tt &a, size_tt &b);    
    void swap(Point<> &a, Point<> &b); 

    // perform a permutation of 'array' of 'len' elements
    // ( we can set 'len' to be any number <= lenght[array])
    // put the first newLen elements of 'array' into newArray
    void rand_Perm(size_tt *array,
                   size_tt len,
                   size_tt *newArray,
                   size_tt newLen);
    
    int min(int a, int b) const { return (a >= b) ? b : a; }

    // convert an int n into a string object.
    std::string itoa(int n);
    // reverse string s in place 
    void reverse(char cstrg[]);
};

#endif

    
