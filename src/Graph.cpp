// Graph.cpp

#include <cmath>
#include <vector>
#include <stdexcept>

#include "Graph.h"


//**************************************************************
//
//      class constructor
//
//**************************************************************
//Graph::Graph()
   
void Graph::clear()
{
    if(adjList){
        for(size_tt v = 0; v < numOfVert + 1; v++)
            delete [] adjList[v];
        delete [] adjList;
        adjList = nullptr;
    }
    if(adjWeight){
        for(size_tt v = 1; v < numOfVert + 1; v++)
            delete [] adjWeight[v];
        delete [] adjWeight;
        adjWeight = nullptr;
    }
    numOfVert = 0;
}

void Graph::from_edge_list(size_tt n,
                           const std::vector<std::pair<size_tt, size_tt>> &edges,
                           const std::vector<coord_t> *weights)
{
    clear();
    numOfVert = n;
    adjList = new size_tt*[numOfVert + 1];
    adjWeight = nullptr;
    size_tt *degs = new size_tt[numOfVert];
    for(size_tt i = 0; i < numOfVert; i++)
        degs[i] = 0;

    bool useWeights = weights && (weights->size() == edges.size());
    if(weights && !useWeights)
        throw std::runtime_error("weights length must match edges length");

    for(const auto &e : edges){
        size_tt u = e.first;
        size_tt v = e.second;
        if(u >= numOfVert || v >= numOfVert)
            throw std::runtime_error("edge list contains vertex out of range");
        if(u == v)
            continue;
        degs[u] += 1;
        degs[v] += 1;
    }

    adjList[0] = degs;
    for(size_tt vert = 1; vert < numOfVert + 1; vert++)
        adjList[vert] = new size_tt[degs[vert-1]];

    if(useWeights){
        adjWeight = new coord_t*[numOfVert + 1];
        adjWeight[0] = nullptr;
        for(size_tt vert = 1; vert < numOfVert + 1; vert++)
            adjWeight[vert] = new coord_t[degs[vert-1]];
    }

    size_tt *count = new size_tt[numOfVert];
    for(size_tt i = 0; i < numOfVert; i++)
        count[i] = 0;

    for(size_t idx = 0; idx < edges.size(); idx++){
        size_tt u = edges[idx].first;
        size_tt v = edges[idx].second;
        if(u == v)
            continue;
        size_tt pos_u = count[u]++;
        size_tt pos_v = count[v]++;
        adjList[u + 1][pos_u] = v;
        adjList[v + 1][pos_v] = u;
        if(useWeights){
            coord_t w = (*weights)[idx];
            adjWeight[u + 1][pos_u] = w;
            adjWeight[v + 1][pos_v] = w;
        }
    }
    delete [] count;
}

void Graph::from_adj_list(size_tt n,
                          const std::vector<std::vector<size_tt>> &adj,
                          const std::vector<std::vector<coord_t>> *weights)
{
    clear();
    numOfVert = n;
    if(adj.size() != numOfVert)
        throw std::runtime_error("adj list size must match n");

    bool useWeights = weights != nullptr;
    if(useWeights && weights->size() != numOfVert)
        throw std::runtime_error("weight list size must match n");

    adjList = new size_tt*[numOfVert + 1];
    size_tt *degs = new size_tt[numOfVert];
    for(size_tt i = 0; i < numOfVert; i++){
        if(useWeights && (*weights)[i].size() != adj[i].size())
            throw std::runtime_error("weight list entries must match adjacency list sizes");
        degs[i] = static_cast<size_tt>(adj[i].size());
    }

    adjList[0] = degs;
    for(size_tt vert = 1; vert < numOfVert + 1; vert++)
        adjList[vert] = new size_tt[degs[vert-1]];

    if(useWeights){
        adjWeight = new coord_t*[numOfVert + 1];
        adjWeight[0] = nullptr;
        for(size_tt vert = 1; vert < numOfVert + 1; vert++)
            adjWeight[vert] = new coord_t[degs[vert-1]];
    } else {
        adjWeight = nullptr;
    }

    for(size_tt v = 0; v < numOfVert; v++){
        for(size_tt j = 0; j < degs[v]; j++){
            size_tt u = adj[v][j];
            if(u >= numOfVert)
                throw std::runtime_error("adj list contains vertex out of range");
            adjList[v + 1][j] = u;
            if(useWeights){
                coord_t w = (*weights)[v][j];
                adjWeight[v + 1][j] = w;
            }
        }
    }
}

void Graph::twistedTorus( size_tt h, size_tt w,
			  size_tt t1, size_tt t2) {

  numOfVert = w * h;
  adjList = new size_tt*[numOfVert+1];
  adjList[0] = new size_tt[numOfVert];

  //setting degrees for all vertices to 4 and allocating memory
  for(size_tt i=0;i<numOfVert;i++) {
    adjList[0][i]=0;
    adjList[i + 1 ] = new size_tt[4]; //adjList[0][i]];
  }

  size_tt a,b,la,lb;
  for (a = 0; a < h; a++) {
    for (b = 0; b < w; b++) {
      // For each vertex match with left neighbor
      la = (a + t1) % h;
      lb = (b + 1) % w;
      addEdge(a+b*h, la+lb*h);
      
      // and bottom neighbor
      la = (a + 1) % h;
      lb = (b + t2) % w;
      addEdge(a+b*h, la+lb*h);
    }
  }
}


//**************************************************************
//
//	torus()
//
//**************************************************************
void Graph::torus( size_tt h, size_tt w )
{

    
    numOfVert = w * h;
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    //setting degrees for all vertices to 4 and allocating memory
    for(size_tt i=0;i<numOfVert;i++) {
      adjList[0][i]=4;
      adjList[i + 1 ] = new size_tt[4]; //adjList[0][i]];
    }

    //
    // Set adjacency lists for all vertices.
    //

    // interior vertices
    for(size_tt i=1; i < w-1; i++)
        for(size_tt v = 1; v < h-1; v++){
            adjList[ (v + i*h) + 1][0] = v + (i-1)*h;
            adjList[ (v + i*h) + 1][1] = v-1 + i*h;
            adjList[ (v + i*h) + 1][2] = v + (i+1)*h;
            adjList[ (v + i*h) + 1][3] = v+1 + i*h;
        }

    // left vertical edge (without corners) 
    for(size_tt v = 1; v < h-1; v++){
        adjList[ v + 1][0] = v - 1;
        adjList[ v + 1][1] = v + h;
        adjList[ v + 1][2] = v + 1;
    }

    // right vertical edge (without corners) 
    for(size_tt v = 1; v < h-1; v++){
        adjList[ (v + (w-1)*h) + 1][0] = v-1 + (w-1)*h;
        adjList[ (v + (w-1)*h) + 1][1] = v   + (w-2)*h;
        adjList[ (v + (w-1)*h) + 1][2] = v+1 + (w-1)*h;
    }

    // top horizontal edge
    for(size_tt i=1; i < w-1; i++){
        adjList[ i*h + 1][0] = (i-1)*h;
        adjList[ i*h + 1][1] = 1 + i*h;
        adjList[ i*h + 1][2] = (i+1)*h;
    }

    // bottom horizontal edge
    for(size_tt i=1; i < w-1; i++){
        adjList[ h + i*h ][0] = h-1 + (i-1)*h;
        adjList[ h + i*h ][1] = h-2 + i*h;
        adjList[ h + i*h ][2] = h-1 + (i+1)*h;
    }

    // FOUR CORNERS

    // top left
    adjList[1][0] = h;
    adjList[1][1] = 1;

    // bottom left
    adjList[h][0] = 2*h-1;
    adjList[h][1] = h-2;

    // top right
    adjList[(w-1)*h + 1][0] = (w-2)*h;
    adjList[(w-1)*h + 1][1] = 1 + (w-1)*h;

    // bottom right
    adjList[ h + (w-1)*h ][0] = h-1 + (w-2)*h;
    adjList[ h + (w-1)*h ][1] = h-2 + (w-1)*h;

    // making the leftmost and rightmost columns adjacent
      for(size_tt j=0;j<h;j++){
	int qqq=(w-1)*h+j;
	adjList[j+1][adjList[0][j]-1]=qqq;
	adjList[qqq+1][adjList[0][j]-1]=j;
        }

      // making the top and bottom row adjacent
      int twist=1;
      for(size_tt j=0;j<numOfVert;j+=h){
          int tw=(h*w+j+h*twist-1)%(h*w);
	  adjList[j+1][3]=tw;
	  adjList[tw+1][3]=j;
      }
}


//**************************************************************
//
//	Method name : pow2
//
//	Description : sequencial doubling of a graph
//
//**************************************************************
void Graph::pow2(size_tt exp)
{

    for(size_tt i = 0; i < exp; i++)
        double_Graph();
    
}
    
//**************************************************************
//
//	Method name : hypercube
//
//	Description : cube of dimension dim
//
//**************************************************************
void Graph::hyper_Cube( size_tt dim )
{

    path_Graph( 2 );
    for(size_tt i = 0; i < dim - 1; i++)
        double_Graph();
    
}

//**************************************************************
//
//	Method name : double_Graph
//
//	Description : create an adjacency list of a graph obtained
//      from the given one G by creating a clone G' of G and joining
//      each vertex of G' with its twin brother in G.
//      
//**************************************************************
void Graph::double_Graph() 
{
    size_tt initSize = numOfVert;
    numOfVert = 2 * initSize;


    size_tt **aList = new size_tt*[numOfVert+1];
    aList[0] = new size_tt[numOfVert];

    // updating degrees
    for(size_tt vert = 0; vert < initSize; vert++)
    {
        aList[ 0 ][ vert ] = adjList[0][ vert ] + 1;
        aList[ 0 ][ vert + initSize ] = aList[0][ vert ];
    }

    // allocating memory and setting adjacency lists of all vertices
    for(size_tt vert = 0; vert < initSize; vert++)
    {// memory allocation
        aList[ vert+1 ] = new size_tt[aList[0][ vert ]];
        aList[ vert+1+initSize ] = new size_tt[aList[0][ vert ]];
            // enlisting a twin brother
        aList[ vert+1 ][0] = initSize + vert;
        aList[ vert+1+initSize ][0] = vert;
        
        for(size_tt adjVert = 1; adjVert < aList[0][ vert ];
            adjVert++)
        {
            aList[ vert+1 ][ adjVert ] =
                adjList[ vert+1 ][ adjVert-1];

            aList[ vert+1+initSize ][ adjVert ] =
                adjList[ vert+1 ][ adjVert-1 ] + initSize;
        }
        
    }
    for(size_tt v = 0; v <= numOfVert; v++)
        delete [] adjList[v];
    delete [] adjList;
    adjList = aList;
    
}

//**************************************************************
//
//	mesh()
//
//**************************************************************
void Graph::mesh( size_tt height )
{

    mesh_Graph(height,height);
    
}

//**************************************************************
//
//	Method name : mesh_Graph
//
//      creating a rectangular w x h mesh
//      
//**************************************************************
void Graph::mesh_Graph(size_tt h, size_tt w)
{

    numOfVert = w * h;
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    // Set degrees of the first column.
    adjList[0][0]=2;
    adjList[0][h-1]=2;
    for(size_tt k=1; k < h-1; k++)
        adjList[0][k]=3;

    
    // memory alloc for the adj list of the elements of the first column
    for(size_tt v = 0; v < h; v++)
        adjList[v+1] = new size_tt[adjList[0][v]];

    // setting degrees and memory alloc of the elements of the remaining columns
    for(size_tt i=1; i<w; i++){
        if( i == w-1 ){
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v];
                adjList[ (v + i*h) + 1 ] = new size_tt[adjList[0][v + i*h]];
            }
        } else {
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v]+1;
                adjList[ v+1+ i*h ] = new size_tt[adjList[0][ v + i*h]];
	    }
	}
    }
    //
    // Set adjacency lists for all vertices.
    //

    // interior vertices
    for(size_tt i=1; i < w-1; i++)
        for(size_tt v = 1; v < h-1; v++){
            adjList[ (v + i*h) + 1][0] = v + (i-1)*h;
            adjList[ (v + i*h) + 1][1] = v-1 + i*h;
            adjList[ (v + i*h) + 1][2] = v + (i+1)*h;
            adjList[ (v + i*h) + 1][3] = v+1 + i*h;
        }

    // left vertical edge (without corners) 
    for(size_tt v = 1; v < h-1; v++){
        adjList[ v + 1][0] = v - 1;
        adjList[ v + 1][1] = v + h;
        adjList[ v + 1][2] = v + 1;
    }

    // right vertical edge (without corners) 
    for(size_tt v = 1; v < h-1; v++){
        adjList[ (v + (w-1)*h) + 1][0] = v-1 + (w-1)*h;
        adjList[ (v + (w-1)*h) + 1][1] = v   + (w-2)*h;
        adjList[ (v + (w-1)*h) + 1][2] = v+1 + (w-1)*h;
    }

    // top horizontal edge
    for(size_tt i=1; i < w-1; i++){
        adjList[ i*h + 1][0] = (i-1)*h;
        adjList[ i*h + 1][1] = 1 + i*h;
        adjList[ i*h + 1][2] = (i+1)*h;
    }

    // bottom horizontal edge
    for(size_tt i=1; i < w-1; i++){
        adjList[ h + i*h ][0] = h-1 + (i-1)*h;
        adjList[ h + i*h ][1] = h-2 + i*h;
        adjList[ h + i*h ][2] = h-1 + (i+1)*h;
    }

    // FOUR CORNERS

    // top left
    adjList[1][0] = h;
    adjList[1][1] = 1;

    // bottom left
    adjList[h][0] = 2*h-1;
    adjList[h][1] = h-2;

    // top right
    adjList[(w-1)*h + 1][0] = (w-2)*h;
    adjList[(w-1)*h + 1][1] = 1 + (w-1)*h;

    // bottom right
    adjList[ h + (w-1)*h ][0] = h-1 + (w-2)*h;
    adjList[ h + (w-1)*h ][1] = h-2 + (w-1)*h;
    
    
}


//**************************************************************
//
//	cylinder_Graph()
//
//**************************************************************
void Graph::square_Cylinder( size_tt h, size_tt w)
{

    numOfVert = w * h;
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    // Set degrees of the first column.
    adjList[0][0]=3;
    adjList[0][h-1]=3;
    for(size_tt k=1; k < h-1; k++)
        adjList[0][k]=4;

    
    // memory alloc for the adj list of the elements of the first column
    for(size_tt v = 0; v < h; v++)
        adjList[v+1] = new size_tt[adjList[0][v]];

    // setting degrees and memory alloc of the elements of the remaining columns
    for(size_tt i=1; i<w; i++){
        if( i == w-1 ){
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v];
                adjList[ (v + i*h) + 1 ] = new size_tt[adjList[0][v + i*h]];
            }
        } else {
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v];
                adjList[ v+1+ i*h ] = new size_tt[adjList[0][ v + i*h]];
	    }
	}
    }
    //
    // Set adjacency lists for all vertices.
    //

    // interior vertices
    for(size_tt i=1; i < w-1; i++)
        for(size_tt v = 1; v < h-1; v++){
            adjList[ (v + i*h) + 1][0] = v + (i-1)*h;
            adjList[ (v + i*h) + 1][1] = v-1 + i*h;
            adjList[ (v + i*h) + 1][2] = v + (i+1)*h;
            adjList[ (v + i*h) + 1][3] = v+1 + i*h;
        }

    // left vertical edge (without corners) 
    for(size_tt v = 1; v < h-1; v++){
        adjList[ v + 1][0] = v - 1;
        adjList[ v + 1][1] = v + h;
        adjList[ v + 1][2] = v + 1;
    }

    // right vertical edge (without corners) 
    for(size_tt v = 1; v < h-1; v++){
        adjList[ (v + (w-1)*h) + 1][0] = v-1 + (w-1)*h;
        adjList[ (v + (w-1)*h) + 1][1] = v   + (w-2)*h;
        adjList[ (v + (w-1)*h) + 1][2] = v+1 + (w-1)*h;
    }

    // top horizontal edge
    for(size_tt i=1; i < w-1; i++){
        adjList[ i*h + 1][0] = (i-1)*h;
        adjList[ i*h + 1][1] = 1 + i*h;
        adjList[ i*h + 1][2] = (i+1)*h;
    }

    // bottom horizontal edge
    for(size_tt i=1; i < w-1; i++){
        adjList[ h + i*h ][0] = h-1 + (i-1)*h;
        adjList[ h + i*h ][1] = h-2 + i*h;
        adjList[ h + i*h ][2] = h-1 + (i+1)*h;
    }

    // FOUR CORNERS

    // top left
    adjList[1][0] = h;
    adjList[1][1] = 1;

    // bottom left
    adjList[h][0] = 2*h-1;
    adjList[h][1] = h-2;

    // top right
    adjList[(w-1)*h + 1][0] = (w-2)*h;
    adjList[(w-1)*h + 1][1] = 1 + (w-1)*h;

    // bottom right
    adjList[ h + (w-1)*h ][0] = h-1 + (w-2)*h;
    adjList[ h + (w-1)*h ][1] = h-2 + (w-1)*h;
    


      for(size_tt j=0;j<h;j++){
	int qqq=(w-1)*h+j;
	adjList[j+1][adjList[0][j]-1]=qqq;
	adjList[qqq+1][adjList[0][j]-1]=j;
        }

}

//**************************************************************
//
//	complete_Graph()
//
//**************************************************************
void Graph::complete_Graph( size_tt _numOfVert )
{
    numOfVert = _numOfVert;

    // the i-th element of the first row holds numbers of adjacent
    // vertices in the i-th row of the adjacency list

    // initializing adjList
    adjList   = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    for(size_tt vert = 0; vert < numOfVert; vert++)
        adjList[ 0 ][ vert ] = numOfVert - 1;

    for(size_tt vert = 0; vert < numOfVert; vert++)
        adjList[vert+1] = new size_tt[numOfVert - 1];
    
    // create an auxiliary array, which holds numbers 0..numOfVert-1
    size_tt *vertices = new size_tt[numOfVert];
    for(size_tt i = 0; i < numOfVert; i++)
        vertices[i] = i;

    size_tt index = 0;
    
    for(size_tt vert = 0; vert < numOfVert; vert++)
    {   // locate vert in 'vertices' and then swap it with the last
        // element of 'vertices'
        index = 0;
        while( vertices[index] !=  vert )
            ++index;
        swap(vertices[index], vertices[numOfVert-1]);
    
        for(size_tt adjVert = 0; adjVert < numOfVert - 1;
            adjVert++)
            adjList[ vert+1 ][ adjVert ] = vertices[adjVert]; 
    }
    delete [] vertices;

}

//**************************************************************
//
//      path_Graph
//
//**************************************************************
void Graph::path_Graph( size_tt _numOfVert )
{

    numOfVert = _numOfVert;

    // initializing adjList
    adjList   = new size_tt*[numOfVert + 1];
    adjList[0] = new size_tt[numOfVert];

    // setting degrees
    // vertex 0 and numOfVert-1 must be dealt with separately
    adjList[0][0] = 1;
    adjList[0][numOfVert-1] = 1;
    for(size_tt vert = 1; vert < numOfVert-1; vert++)
    adjList[ 0 ][ vert ] = 2;

    // setting adjLists for each vertex
    // vertex 0 and numOfVert-1 must be considered separately
    adjList[1]    = new size_tt[1];
    adjList[1][0] = 1;
    
    adjList[numOfVert]    = new size_tt[1];
    adjList[numOfVert][0] = numOfVert-2;
        
    for(size_tt vert = 1; vert < numOfVert-1; vert++){
        adjList[ vert+1 ] = new size_tt[adjList[ 0 ][ vert ]];
        for(size_tt adjVert = 0; adjVert < adjList[ 0 ][ vert ];
            adjVert++) {
            adjList[ vert+1 ][ adjVert ] = vert - 1 + 2*adjVert; 
        }
    }
}

//**************************************************************
//
//	cycle_Graph()
//
//**************************************************************
void Graph::cycle_Graph( size_tt _numOfVert )
{

    numOfVert = _numOfVert;

    // initializing adjList
    adjList   =  new size_tt*[numOfVert + 1];
    adjList[0] = new size_tt[numOfVert];

    // setting degrees
    for(size_tt vert = 0; vert < numOfVert; vert++)
    adjList[ 0 ][ vert ] = 2;


    // allocating memory
    for(size_tt vert = 0; vert < numOfVert; vert++)
        adjList[ vert+1 ] = new size_tt[adjList[ 0 ][ vert ]];

    // setting adjacency lists for each vertex
    // vertices 0 and numOfVert-1 must be dealt with separately
    adjList[ 1 ][ 0 ] = numOfVert - 1;
    adjList[ 1 ][ 1 ] = 1;

    adjList[numOfVert][0] = numOfVert-2;
    adjList[numOfVert][1] = 0;

    
    for(size_tt vert = 1; vert < numOfVert-1; vert++)
        for(size_tt adjVert = 0; adjVert < adjList[ 0 ][ vert ];
            adjVert++){
            adjList[ vert+1 ][ adjVert ] = vert-1 + 2*adjVert; 
        }

}

//****************************************************************
//
//	rand_cpt_Graph()
//
//****************************************************************
size_tt ** Graph::rand_cpt_Graph( size_tt _numOfVert )
{
  sfast_Rand(static_cast<unsigned int>(time(0)));

  numOfVert = _numOfVert;

  std::vector<size_tt> adj(numOfVert-1); // vector whose entries count the number of elements in the corresponding rows of adjList

  for(size_tt index = 0; index < numOfVert-1; index++){
    // for each vertex vert between 0 and numOfVert-2 generate
    // a random number between 1 and numOfVert - 1 - vert.
    adj[index] = 1 + fast_Rand() % (numOfVert - index - 1);
  }

  
  // initializing adjList
  size_tt ** cptAdjList = new size_tt*[numOfVert];
  cptAdjList[0] = new size_tt[numOfVert - 1];

  for(size_tt adjVertex = 0; adjVertex < numOfVert - 1; adjVertex++)
    cptAdjList[0][adjVertex] = adj[adjVertex];
    
  size_tt *array = new size_tt[numOfVert];

  for(size_tt vertex = 1; vertex < numOfVert; vertex++) {
    // generate adj[ vertex ] random size_tt without repetition
    // in the range vertex .. numOfVert-1
    for(size_tt i = 0; i < numOfVert - vertex; i++)
      array[i] = i + vertex;
        
    rand_Perm(array, (size_tt)(numOfVert - vertex),
              cptAdjList[vertex], adj[vertex - 1]);
  }

  delete [] array;

  return cptAdjList;
}    

//****************************************************************
//
//	rand_Graph()
//
//****************************************************************
void Graph::rand_Graph( size_tt h, size_tt w)
{

  sfast_Rand(static_cast<unsigned int>(time(0)));
  numOfVert=h;
  size_tt type=0;
  if(fast_Rand() % 2 == 1)
    type=1;
  size_tt temp=0;
  size_tt total=0;
  while(total<=h) {
    total=total+(size_tt)pow(2,temp);
    temp++;
  }
  if (total> h)
    total=total-(size_tt)pow(2,temp-1);
  size_tt first_leaf=total-(size_tt)pow(2,temp-2);
  numOfVert=total;

  
  if (type==0){
    size_tt base=2;
    adjList=new size_tt*[numOfVert+1];
    adjList[0]=new size_tt[numOfVert];
    adjList[0][0]=base;
    adjList[1] = new size_tt[h-1];
    
    for(size_tt i=1;i<first_leaf;i++) {
      adjList[0][i]=base+1;
      adjList[i+1] = new size_tt[h-1];
    }
    for(size_tt i=first_leaf;i<numOfVert;i++) {
      adjList[0][i]=1;
      adjList[i+1] = new size_tt[h-1];
    }
    
    //setting adjacencies for the root  
    for (size_tt i=0;i<base;i++)
      adjList[1][i]=i+1;
    
    //setting adjacencies for the internal nodes
    
    for(size_tt i=1;i<first_leaf;i++){
      adjList[i+1][0]=(size_tt)((i-1)/base);
      for(size_tt j=1;j<=base;j++)
	adjList[i+1][j]=base*i+j;
    }
    
    //setting adjacencies for the leafs
    for(size_tt i=first_leaf;i<numOfVert;i++)
      adjList[i+1][0]=(size_tt)((i-1)/base);
    
  }
  
  else{
    
    // initializing adjList
    adjList   = new size_tt*[numOfVert + 1];
    adjList[0] = new size_tt[numOfVert];
    
    // setting degrees
    // vertex 0 and numOfVert-1 must be dealt with separately
    adjList[0][0] = 1;
    adjList[0][numOfVert-1] = 1;
    for(size_tt vert = 1; vert < numOfVert-1; vert++)
      adjList[ 0 ][ vert ] = 2;
    
    // setting adjLists for each vertex
    // vertex 0 and numOfVert-1 must be considered separately
    adjList[1]    = new size_tt[1];
    adjList[1][0] = 1;
    
    adjList[numOfVert]    = new size_tt[1];
    adjList[numOfVert][0] = numOfVert-2;
    
    for(size_tt vert = 1; vert < numOfVert-1; vert++){
      adjList[ vert+1 ] = new size_tt[numOfVert-1];
      adjList[vert+1][0]=vert-1;
      adjList[vert+1][1]=vert+1;
    } 
  }

  for(size_tt i=0;i<numOfVert;i++){
    for(size_tt j=0;j<i;j++){
      if ((i!=j)&&(i!=(j-1))&&(i!=j+1)) {
	size_tt th = fast_Rand() % (numOfVert * numOfVert);
	if ((i<5)&&(th<=(w*w))) {
	  adjList[i+1][adjList[0][i]++]=j;
	  adjList[j+1][adjList[0][j]++]=i;
	}
	else if ((i>h-4)&&(j>h-4)&&(th<=(w*w))) {
	  adjList[i+1][adjList[0][i]++]=j;
	  adjList[j+1][adjList[0][j]++]=i;
	}
	else if (th<=w) {
	  adjList[i+1][adjList[0][i]++]=j;
	  adjList[j+1][adjList[0][j]++]=i;
	}
      }
    }
  }
}

//****************************************************************
//
//	Method name : rand_Perm
//
//	Description : perform a permutation of 'array' of 'len'
//      elements (we can set 'len' to be any number <= length(array))
//      and put the first newLen elements of 'array' into newArray.
//
//****************************************************************
 void Graph::rand_Perm(size_tt *array,
                      size_tt len,
                      size_tt *newArray,
                      size_tt newLen)
{
    size_tt l = (newLen == 0) ? len : newLen;
    
    for( size_tt i = 0; i < l; i++)
        swap(array[i], array[i + (fast_Rand() % (len - i))]);
    
    for( size_tt i = 0; i < newLen; i++)
        newArray[i] = array[i];
}

//****************************************************************
//
//	swap()
//
//	swap two int, coord_t or Point<>
//
//****************************************************************
void Graph::swap(int &a, int &b)
{
    int temp = a;
    a = b;
    b = temp;
}

void Graph::swap(size_tt &a, size_tt &b)
{
    size_tt temp = a;
    a = b;
    b = temp;
}

void Graph::swap(Point<> &a, Point<> &b)
{
    Point<> temp = a;
    a = b;
    b = temp;
}
//****************************************************************
//
//	Method name : fast_Rand
//
//      fast random number generator from "Numerical Recipes in C"
//      section "An even quicker generator" p. 284.  
//
//****************************************************************
unsigned long idum;
void Graph::sfast_Rand(unsigned int seed) const { idum = seed; }
unsigned long Graph::fast_Rand() const {
    idum = 1664525L*idum + 1013904223L;
    return idum;
}

//****************************************************************
//
//    get_Max_Deg_Vert()
//
//    minimal degree of the graph
//
//****************************************************************
size_tt Graph::get_Max_Deg_Vert(size_tt &max_deg, size_tt &min_deg) const
{
    size_tt max_deg_vert = 0;
    max_deg = 0;
    min_deg = numOfVert;
    
    for(size_tt vert = 0; vert < numOfVert; vert++){
        if( max_deg < adjList[0][vert] ){
            max_deg = adjList[0][vert];
            max_deg_vert = vert;
        }
        if( min_deg > adjList[0][vert] ){
            min_deg = adjList[0][vert];
        }
    }
    return max_deg_vert;
}

//**************************************************************
//
//	square_Torus()
//
//**************************************************************
void Graph::square_Torus( size_tt base )
{

    cycle_Graph( base );
    for(size_tt i = 0; i < 2; i++)
        double_Graph();
    
}

//**************************************************************
//
//      tree()
//
//**************************************************************
void Graph::tree(size_tt depth,size_tt base)
{
  //make sure that the tree is at least binary at depth 2

  if (depth<2) 
    depth=2;
  if (base<2) 
    base=2;
  
  //compute the total number of vertices of the tree
  size_tt total=0;
  for(size_tt i=0;i<depth;i++)
    total = total + (size_tt)pow(base,i);
  numOfVert = total;

  adjList=new size_tt*[numOfVert+1];
  adjList[0]=new size_tt[numOfVert];
  adjList[0][0]=base;
  adjList[1] = new size_tt[base];
  
  size_tt first_leaf=numOfVert-(size_tt)pow(base,depth-1);

  for(size_tt i=1;i<first_leaf;i++) {
    adjList[0][i]=base+1;
    adjList[i+1] = new size_tt[base+1];
  }
  for(size_tt i=first_leaf;i<numOfVert;i++) {
    adjList[0][i]=1;
    adjList[i+1] = new size_tt[1];
  }

  //setting adjacencies for the root  
  for (size_tt i=0;i<base;i++)
    adjList[1][i]=i+1;

  //setting adjacencies for the internal nodes

  for(size_tt i=1;i<first_leaf;i++){
    adjList[i+1][0]=(size_tt)((i-1)/base);
    for(size_tt j=1;j<=base;j++)
      adjList[i+1][j]=base*i+j;
  }
		 
  //setting adjacencies for the leafs
  for(size_tt i=first_leaf;i<numOfVert;i++)
    adjList[i+1][0]=(size_tt)((i-1)/base);


}

//**************************************************************
//
//      meshX()
//
//**************************************************************
void Graph::meshX( size_tt height )
{
  size_tt w=height;
  size_tt h=height;
    numOfVert = w * h;
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    // Set degrees of the first column.
    adjList[0][0]=2;
    adjList[0][h-1]=2;
    for(size_tt k=1; k < h-1; k++)
        adjList[0][k]=3;

    
    // memory alloc for the adj list of the elements of the first column
    for(size_tt v = 0; v < h; v++)
        adjList[v+1] = new size_tt[adjList[0][v]];

    // setting degrees and memory alloc of the elements of the remaining columns
    for(size_tt i=1; i<w; i++){
        if( i == w-1 ){
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v];
                adjList[ (v + i*h) + 1 ] = new size_tt[adjList[0][v + i*h]];
            }
        } else {
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v]+1;
                adjList[ v+1+ i*h ] = new size_tt[adjList[0][ v + i*h]];
	    }
	}
    }
    //
    // Set adjacency lists for all vertices.
    //

    // interior vertices
    for(size_tt i=1; i < w-1; i++)
        for(size_tt v = 1; v < h-1; v++){
            adjList[ (v + i*h) + 1][0] = v + (i-1)*h;
            adjList[ (v + i*h) + 1][1] = v-1 + i*h;
            adjList[ (v + i*h) + 1][2] = v + (i+1)*h;
            adjList[ (v + i*h) + 1][3] = v+1 + i*h;
        }

    // left vertical edge (without corners) 
    for(size_tt v = 1; v < h-1; v++){
        adjList[ v + 1][0] = v - 1;
        adjList[ v + 1][1] = v + h;
        adjList[ v + 1][2] = v + 1;
    }

    // right vertical edge (without corners) 
    for(size_tt v = 1; v < h-1; v++){
        adjList[ (v + (w-1)*h) + 1][0] = v-1 + (w-1)*h;
        adjList[ (v + (w-1)*h) + 1][1] = v   + (w-2)*h;
        adjList[ (v + (w-1)*h) + 1][2] = v+1 + (w-1)*h;
    }

    // top horizontal edge
    for(size_tt i=1; i < w-1; i++){
        adjList[ i*h + 1][0] = (i-1)*h;
        adjList[ i*h + 1][1] = 1 + i*h;
        adjList[ i*h + 1][2] = (i+1)*h;
    }

    // bottom horizontal edge
    for(size_tt i=1; i < w-1; i++){
        adjList[ h + i*h ][0] = h-1 + (i-1)*h;
        adjList[ h + i*h ][1] = h-2 + i*h;
        adjList[ h + i*h ][2] = h-1 + (i+1)*h;
    }

    // FOUR CORNERS

    // top left
    adjList[1][0] = h;
    adjList[1][1] = 1;

    // bottom left
    adjList[h][0] = 2*h-1;
    adjList[h][1] = h-2;

    // top right
    adjList[(w-1)*h + 1][0] = (w-2)*h;
    adjList[(w-1)*h + 1][1] = 1 + (w-1)*h;

    // bottom right
    adjList[ h + (w-1)*h ][0] = h-1 + (w-2)*h;
    adjList[ h + (w-1)*h ][1] = h-2 + (w-1)*h;
    
    

    adjList[0][0]=5;
    adjList[0][h-1]=5;
    adjList[0][h*w-1]=5;
    adjList[0][h*w-w]=5;
    
    adjList[1][2]=h-1;
    adjList[1][3]=h*w-1;
    adjList[1][4]=h*w-w;
    
    adjList[h][2]=0;
    adjList[h][3]=h*w-1;
    adjList[h][4]=h*w-w;

    adjList[h*w][2]=h-1;
    adjList[h*w][3]=0;
    adjList[h*w][4]=h*w-w;

    adjList[h*w-w+1][2]=h-1;
    adjList[h*w-w+1][3]=h*w-1;
    adjList[h*w-w+1][4]=0;


    

}

//**************************************************************
//
//      meshT()
//
//**************************************************************
void Graph::meshT(size_tt h)
{
  numOfVert=(h+1)*h/2;
  size_tt index=0;
  adjList = new size_tt*[numOfVert+1];
  adjList[0]=new size_tt[numOfVert];
  adjList[1]=new size_tt[2];
  adjList[0][0]=2;
  adjList[1][0]=1;
  adjList[1][1]=2;

  for(size_tt i=2;i<=h;i++) {
    for(size_tt j=1;j<=i;j++){
      index=index+1;
      if (i!=h){
	if (j==1) {
	  adjList[0][index]=4;
	  adjList[index+1]=new size_tt[4];
	  adjList[index+1][0]=index-(i-1);
	  adjList[index+1][1]=index+1;
	  adjList[index+1][2]=index+i;
	  adjList[index+1][3]=index+i+1;
	}
	else if(j==i) {
	  adjList[0][index]=4;
	  adjList[index+1]=new size_tt[4];
	  adjList[index+1][0]=index-i;
	  adjList[index+1][1]=index-1;
	  adjList[index+1][2]=index+i;
	  adjList[index+1][3]=index+i+1;	
	}
	else {
	  adjList[0][index]=6;
	  adjList[index+1]=new size_tt[6];
	  adjList[index+1][0]=index-i;
	  adjList[index+1][1]=index-i+1;
	  adjList[index+1][2]=index-1;
	  adjList[index+1][3]=index+1;
	  adjList[index+1][4]=index+i;
	  adjList[index+1][5]=index+i+1;
	}
      }
      else {
	if(j==1) {	
	  adjList[0][index]=2;
	  adjList[index+1]=new size_tt[2];
	  adjList[index+1][0]=index-(i-1);
	  adjList[index+1][1]=index+1;
	}
	else if(j==i) {
	  adjList[0][index]=2;
	  adjList[index+1]=new size_tt[2];
	  adjList[index+1][0]=index-i;
	  adjList[index+1][1]=index-1;
	}
	else {
	  adjList[0][index]=4;
	  adjList[index+1]=new size_tt[4];
	  adjList[index+1][0]=index-i;
	  adjList[index+1][1]=index-i+1;
	  adjList[index+1][2]=index-1;
	  adjList[index+1][3]=index+1;
	}
      }
    }
  }
}

//**************************************************************
//
//      meshTX()
//
//**************************************************************
void Graph::meshTX(size_tt h)
{
  numOfVert=(h+1)*h/2+1;
  size_tt index=0;

  adjList = new size_tt*[numOfVert+1];
  adjList[0]=new size_tt[numOfVert];
  adjList[0][0]=3;

  adjList[1]=new size_tt[3];
  adjList[1][0]=1;
  adjList[1][1]=2;
  adjList[1][2]=numOfVert-1;
  adjList[0][numOfVert-1]=3;
  adjList[numOfVert]=new size_tt[3];
  adjList[numOfVert][0]=0;
  adjList[numOfVert][1]=numOfVert-h-1;
  adjList[numOfVert][2]=numOfVert-2;

  for(size_tt i=2;i<=h;i++) {
    for(size_tt j=1;j<=i;j++){
      index=index+1;
      if (i!=h){
	if (j==1) {
	  adjList[0][index]=4;
	  adjList[index+1]=new size_tt[4];
	  adjList[index+1][0]=index-(i-1);
	  adjList[index+1][1]=index+1;
	  adjList[index+1][2]=index+i;
	  adjList[index+1][3]=index+i+1;
	}
	else if(j==i) {
	  adjList[0][index]=4;
	  adjList[index+1]=new size_tt[4];
	  adjList[index+1][0]=index-i;
	  adjList[index+1][1]=index-1;
	  adjList[index+1][2]=index+i;
	  adjList[index+1][3]=index+i+1;	
	}
	else {
	  adjList[0][index]=6;
	  adjList[index+1]=new size_tt[6];
	  adjList[index+1][0]=index-i;
	  adjList[index+1][1]=index-i+1;
	  adjList[index+1][2]=index-1;
	  adjList[index+1][3]=index+1;
	  adjList[index+1][4]=index+i;
	  adjList[index+1][5]=index+i+1;
	}
      }
      else {
	if(j==1) {	
	  adjList[0][index]=3;
	  adjList[index+1]=new size_tt[3];
	  adjList[index+1][0]=index-(i-1);
	  adjList[index+1][1]=index+1;
	  adjList[index+1][2]=numOfVert-1;
	}
	else if(j==i) {
	  adjList[0][index]=3;
	  adjList[index+1]=new size_tt[3];
	  adjList[index+1][0]=index-i;
	  adjList[index+1][1]=index-1;
	  adjList[index+1][2]=numOfVert-1;
	}
	else {
	  adjList[0][index]=4;
	  adjList[index+1]=new size_tt[4];
	  adjList[index+1][0]=index-i;
	  adjList[index+1][1]=index-i+1;
	  adjList[index+1][2]=index-1;
	  adjList[index+1][3]=index+1;
	}
      }
    }
  }
}

inline void Graph::addEdge(size_tt a, size_tt b)
{
  adjList[a+1][adjList[0][a]] = b;
  adjList[b+1][adjList[0][b]] = a;
  adjList[0][a]++;
  adjList[0][b]++;
}

//**************************************************************
//
//      sierpinski_recurse3D()
//
//**************************************************************
static size_tt vert = 4;
size_tt Graph::sierpinski_recurse3D(int maxLevel,
				   int currentLevel,
				   size_tt a, size_tt b, 
				   size_tt c,
				   size_tt d,
				   bool flag)
{ 
  if (currentLevel >= maxLevel) {
    if (flag) {
      addEdge(a,b);
      addEdge(a,c);
      addEdge(a,d);
      addEdge(b,c);
      addEdge(b,d);
      addEdge(c,d);
    }
  } else {
    size_tt e = vert++;
    size_tt f = vert++;
    size_tt g = vert++;
    size_tt h = vert++;
    size_tt i = vert++;
    size_tt j = vert++;
    sierpinski_recurse3D(maxLevel,currentLevel+1,a,g,e,f,flag);
    sierpinski_recurse3D(maxLevel,currentLevel+1,e,b,i,h,flag);
    sierpinski_recurse3D(maxLevel,currentLevel+1,f,c,j,h,flag);
    sierpinski_recurse3D(maxLevel,currentLevel+1,g,d,i,j,flag);
  }
  return vert;
}
 
void Graph::sierpinski_recurse(int maxLevel,
			      int currentLevel,
			      size_tt a, size_tt b, size_tt c)
{
  static size_tt vert = 3;
  if (currentLevel >= maxLevel) {
    addEdge(a,b);
    addEdge(a,c);
    addEdge(b,c);
  } else {
    size_tt d = vert++;
    size_tt e = vert++;
    size_tt f = vert++;
    sierpinski_recurse(maxLevel,currentLevel+1,a,d,e);
    sierpinski_recurse(maxLevel,currentLevel+1,d,b,f);
    sierpinski_recurse(maxLevel,currentLevel+1,e,f,c);
  }
}

//**************************************************************
//
//      sierpinski()
//
//**************************************************************
void Graph::sierpinski(size_tt h, size_tt density)
{
  if (density == 0) {
    size_tt temp=3;
    for(size_tt i=1;i<=h;i++)
      temp=temp+(size_tt)pow(3,i);
    numOfVert=temp;
    adjList = new size_tt*[numOfVert+1];
    adjList[0]=new size_tt[numOfVert];
    for(size_tt i=0;i<numOfVert;i++) {
      adjList[i+1] = new size_tt[4];
      adjList[0][i]=0;
    }
    
    sierpinski_recurse(h,0,0,1,2);
  }
  else {
    numOfVert = sierpinski_recurse3D(h,0,0,1,2,3,false);
    adjList = new size_tt*[numOfVert+1];
    adjList[0]=new size_tt[numOfVert];
    for(size_tt i=0;i<numOfVert;i++) {
      adjList[i+1] = new size_tt[6];
      adjList[0][i]=0;
    }
    vert=4;
    sierpinski_recurse3D(h,0,0,1,2,3, true);
  }
}
