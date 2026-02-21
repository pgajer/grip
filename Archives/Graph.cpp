// Graph.cpp

#include "Graph.hpp"
#include <cmath>
#define DEBUG 0
#define DEBUG_RS 0

//**************************************************************
//
//      class constructor
//
//**************************************************************
//Graph::Graph()
   
void Graph::twistedTorus( size_tt h, size_tt w,
			  size_tt t1, size_tt t2) {
  debug("h = " << h);
  debug("w = " << w);
  debug("t1 = " << t1);
  debug("t2 = " << t2);

  numOfVert = w * h;
  adjList = new size_tt*[numOfVert+1];
  adjList[0] = new size_tt[numOfVert];

  //setting degrees for all vertices to 4 and allocating memory
  for(int i=0;i<numOfVert;i++) {
    adjList[0][i]=0;
    adjList[i + 1 ] = new size_tt[4]; //adjList[0][i]];
  }

  size_tt a,b,la,lb;
  for (a = 0; a < h; a++) {
    for (b = 0; b < w; b++) {
      // For each vertex match with left neighbor
      la = (a + t1) % h;
      lb = (b + 1) % w;
      debug("Adding edge " << a+b*h << " to " << la+lb*h);
      addEdge(a+b*h, la+lb*h);
      
      // and bottom neighbor
      la = (a + 1) % h;
      lb = (b + t2) % w;
      debug("Adding edge " << a+b*h << " to " << la+lb*h);
      addEdge(a+b*h, la+lb*h);
    }
  }
#ifdef __SKIP__
  debug("adjList: ");
  for(size_tt v=0; v< numOfVert; v++){
    cout << v << ": (adjList[0]["<<v<<"]="<<adjList[0][v]<<") : ";
    for(size_tt j=0; j < adjList[0][v]; j++){
      cout << adjList[v+1][j] << ' ';
    }
    cout << endl;
  }
  exit(1);
  debug("Leaving twisted torus");
#endif
}


//**************************************************************
//
//	torus()
//
//**************************************************************
void Graph::torus( size_tt h, size_tt w )
{
#if DEBUG
    debug("Entering torus");
#endif 

    
    numOfVert = w * h;
#if DEBUG_MISH_GRAPH
    debug("h=" << h<<", w="<< w<<", numOfVert=" << numOfVert );
#endif
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    //setting degrees for all vertices to 4 and allocating memory
    for(int i=0;i<numOfVert;i++) {
      adjList[0][i]=4;
      adjList[i + 1 ] = new size_tt[4]; //adjList[0][i]];
    }

    //
    //  SETTING ADJACENCY LISTS OF ALL VICES
    //

    // interior vertices
    for(int i=1; i < w-1; i++)
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
    for(int i=1; i < w-1; i++){
        adjList[ i*h + 1][0] = (i-1)*h;
        adjList[ i*h + 1][1] = 1 + i*h;
        adjList[ i*h + 1][2] = (i+1)*h;
    }

    // bottom horizontal edge
    for(int i=1; i < w-1; i++){
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
      for(int j=0;j<h;j++){
	int qqq=(w-1)*h+j;
	adjList[j+1][adjList[0][j]-1]=qqq;
	adjList[qqq+1][adjList[0][j]-1]=j;
        }

      // making the top and bottom row adjacent
      int twist=1;
      for(int j=0;j<numOfVert;j+=h){
	//	for(int k=0;k<adjList[0][j];k++) {
	//adjList[j+1][k]=adjList[j+1][k];
	//}
	//	if((j%h)==0) {
          int tw=(h*w+j+h*twist-1)%(h*w);
	  adjList[j+1][3]=tw;
	  adjList[tw+1][3]=j;
	  // }
      }
#if DEBUG
    debug("Leaving torus");
#endif 
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
#if DEBUG
    debug("Entering pow2");
#endif 

    for(size_tt i = 0; i < exp; i++)
        double_Graph();
    
#if DEBUG
    debug("Leaving pow2");
#endif 
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
#if DEBUG
    debug("Entering hypercube");
#endif 

    path_Graph( 2 );
    for(size_tt i = 0; i < dim - 1; i++)
        double_Graph();
    
#if DEBUG
    debug("Leaving hypercube");
#endif 
//    return *this;
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
//#define DEBUG
#if DEBUG
    debug("Entering double_Graph");
#endif
    size_tt initSize = numOfVert;
    numOfVert = 2 * initSize;

//    debug("BEFORE\n" << adjList);

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
#if DEBUG     
        debug("vert = " << vert );
#endif
            // enlisting a twin brother
        aList[ vert+1 ][0] = initSize + vert;
        aList[ vert+1+initSize ][0] = vert;
        
        for(size_tt adjVert = 1; adjVert < aList[0][ vert ];
            adjVert++)
        {
#if DEBUG   
            debug("1st, adjVert = " << adjVert );
#endif    
            aList[ vert+1 ][ adjVert ] =
                adjList[ vert+1 ][ adjVert-1];

#if DEBUG     
            debug("2nd, adjVert = " << adjVert );
#endif      
            aList[ vert+1+initSize ][ adjVert ] =
                adjList[ vert+1 ][ adjVert-1 ] + initSize;
        }
        
    }
    for(size_tt v = 0; v <= numOfVert; v++)
        delete [] adjList[v];
    delete [] adjList;
    adjList = aList;
//	debug("AFTER\n" << adjList);
    
#if DEBUG
    debug("*this = " << *this );
    debug("Leaving double_Graph");
#endif
//    return *this;
}

//**************************************************************
//
//	mesh()
//
//**************************************************************
void Graph::mesh( size_tt height )
{
#if DEBUG 
    debug("Entering");
#endif 

    mesh_Graph(height,height);
    
#if DEBUG 
    debug("Leaving");
#endif 
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
#define DEBUG_MISH_GRAPH 0
#if DEBUG_MISH_GRAPH
    debug("Entering mesh_Graph");
#endif

    numOfVert = w * h;
#if DEBUG_MISH_GRAPH
    debug("h=" << h<<", w="<< w<<", numOfVert=" << numOfVert );
#endif
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    //setting degrees of the first colunm
    adjList[0][0]=2;
    adjList[0][h-1]=2;
    for(int k=1; k < h-1; k++)
        adjList[0][k]=3;

//      for(int v=0; v < h; v++)
//          debug("adjList[0]["<<v<<"]="<<adjList[0][v]);
    
    // memory alloc for the adj list of the elements of the first column
    for(size_tt v = 0; v < h; v++)
        adjList[v+1] = new size_tt[adjList[0][v]];

    // setting degrees and memory alloc of the elements of the remaining columns
    for(int i=1; i<w; i++){     
//        debug("i=" << i);
        if( i == w-1 ){
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v];
                adjList[ (v + i*h) + 1 ] = new size_tt[adjList[0][v + i*h]];
//                debug("adjList[0]["<<v+i*h<<"]="<<adjList[0][v + i*h]);
            }
        } else {
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v]+1;
                adjList[ v+1+ i*h ] = new size_tt[adjList[0][ v + i*h]];
//                debug("adjList[0]["<<v+i*h<<"]="<<adjList[0][ v + i*h ]); 
	    }
	}
    }
    //
    //  SETTING ADJACENCY LISTS OF ALL VICES
    //

    // interior vertices
    for(int i=1; i < w-1; i++)
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
    for(int i=1; i < w-1; i++){
        adjList[ i*h + 1][0] = (i-1)*h;
        adjList[ i*h + 1][1] = 1 + i*h;
        adjList[ i*h + 1][2] = (i+1)*h;
    }

    // bottom horizontal edge
    for(int i=1; i < w-1; i++){
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
    
//      debug("adjList: ");
//      for(size_tt v=0; v< numOfVert; v++){
//          cout << v << ": (adjList[0]["<<v<<"]="<<adjList[0][v]<<") : ";
//          for(size_tt j=0; j < adjList[0][v]; j++){
//              cout << adjList[v+1][j] << ' ';
//          }
//          cout << endl;
//      }
//      exit(1);
    
#if DEBUG_MISH_GRAPH 
    debug("Leaving mesh_Graph");
#endif
}


//**************************************************************
//
//	cylinder_Graph()
//
//**************************************************************
void Graph::square_Cylinder( size_tt h, size_tt w)
{
#define DEBUG_MISH_GRAPH 0
#if DEBUG_MISH_GRAPH
    debug("Entering mesh_Graph");
#endif

    numOfVert = w * h;
#if DEBUG_MISH_GRAPH
    debug("h=" << h<<", w="<< w<<", numOfVert=" << numOfVert );
#endif
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    //setting degrees of the first colunm
    adjList[0][0]=3;
    adjList[0][h-1]=3;
    for(int k=1; k < h-1; k++)
        adjList[0][k]=4;

//      for(int v=0; v < h; v++)
//          debug("adjList[0]["<<v<<"]="<<adjList[0][v]);
    
    // memory alloc for the adj list of the elements of the first column
    for(size_tt v = 0; v < h; v++)
        adjList[v+1] = new size_tt[adjList[0][v]];

    // setting degrees and memory alloc of the elements of the remaining columns
    for(int i=1; i<w; i++){     
//        debug("i=" << i);
        if( i == w-1 ){
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v];
                adjList[ (v + i*h) + 1 ] = new size_tt[adjList[0][v + i*h]];
//                debug("adjList[0]["<<v+i*h<<"]="<<adjList[0][v + i*h]);
            }
        } else {
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v];
                adjList[ v+1+ i*h ] = new size_tt[adjList[0][ v + i*h]];
//                debug("adjList[0]["<<v+i*h<<"]="<<adjList[0][ v + i*h ]); 
	    }
	}
    }
    //
    //  SETTING ADJACENCY LISTS OF ALL VICES
    //

    // interior vertices
    for(int i=1; i < w-1; i++)
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
    for(int i=1; i < w-1; i++){
        adjList[ i*h + 1][0] = (i-1)*h;
        adjList[ i*h + 1][1] = 1 + i*h;
        adjList[ i*h + 1][2] = (i+1)*h;
    }

    // bottom horizontal edge
    for(int i=1; i < w-1; i++){
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
    
//      debug("adjList: ");
//      for(size_tt v=0; v< numOfVert; v++){
//          cout << v << ": (adjList[0]["<<v<<"]="<<adjList[0][v]<<") : ";
//          for(size_tt j=0; j < adjList[0][v]; j++){
//              cout << adjList[v+1][j] << ' ';
//          }
//          cout << endl;
//      }
//      exit(1);


      for(int j=0;j<h;j++){
	int qqq=(w-1)*h+j;
	adjList[j+1][adjList[0][j]-1]=qqq;
	adjList[qqq+1][adjList[0][j]-1]=j;
        }

#if DEBUG_MISH_GRAPH 
    debug("Leaving mesh_Graph");
#endif
}

//**************************************************************
//
//	complete_Graph()
//
//**************************************************************
void Graph::complete_Graph( size_tt _numOfVert )
{
#define DEBUG_COMPL 0
#if DEBUG || DEBUG_COMPL
    debug("Entering complete_Graph");
#endif
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

#if DEBUG
    debug("Exiting complete_Graph");
#endif
//    return *this;
}

//**************************************************************
//
//      path_Graph
//
//**************************************************************
void Graph::path_Graph( size_tt _numOfVert )
{
#define DEBUG_PATH 0
#if DEBUG || DEBUG_PATH
    debug("Entering path");
#endif

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
#if DEBUG_PATH
    debug("Exiting path");
#endif
//    return *this;
}

//**************************************************************
//
//	cycle_Graph()
//
//**************************************************************
void Graph::cycle_Graph( size_tt _numOfVert )
{
#define DEBUGa 0
#if DEBUGa
    debug("Entering cycle");
#endif

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

#if DEBUGa
    debug("Exiting cycle");
#endif
//    return *this;
}

//****************************************************************
//
//	rand_cpt_Graph()
//
//****************************************************************
size_tt ** Graph::rand_cpt_Graph( size_tt _numOfVert )
{
#if DEBUG
    debug("Entering rand_cpt_Graph" );
#endif

    srand(time(0));

    numOfVert = _numOfVert;
    
    size_tt adj[numOfVert-1]; // array whose entries count the number of
                          // elements in the corresponding rows of
                          // adjList
    
    for(size_tt index = 0; index < numOfVert-1; index++){
        // for each vertex vert between 0 and numOfVert-2 generate
        // a random number between 1 and numOfVert - 1 - vert.
        adj[ index ] = 1 + rand() % (numOfVert - index - 1);
        
#if DEBUG
        debug("Inside for loop, with index = "<<
            index << endl <<"adj[" << index << "] = "
             << adj[index] );
#endif
    }

  
    // initializing adjList
    size_tt ** cptAdjList = new size_tt*[numOfVert];
    cptAdjList[0] = new size_tt[ numOfVert - 1];

    for(size_tt adjVertex = 0; adjVertex < numOfVert - 1; adjVertex++)
        cptAdjList[ 0 ][ adjVertex ] = adj[ adjVertex ];
    

    for(size_tt vertex = 1; vertex < numOfVert; vertex++) {
        // generate adj[ vertex ] random size_tt without repetition
        // in the range vertex .. numOfVert-1
        size_tt *array = new size_tt[numOfVert - vertex];
        for(size_tt i = 0; i < numOfVert - vertex; i++)
            array[i] = i + vertex;
        
        rand_Perm(array, (size_tt)(numOfVert - vertex),
                  cptAdjList[ vertex ], adj[ vertex - 1]);

        delete [] array;
    }    
#if DEBUG
    debug("Exiting rand_cpt_Graph" );
#endif
        return cptAdjList;
}    

//****************************************************************
//
//	rand_Graph()
//
//****************************************************************
void Graph::rand_Graph( size_tt h, size_tt w)
{

  srand(time(0));
  numOfVert=h;
  size_tt type=0;
  if(rand()%2==1)
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
  debug("total="<<total<<", numOfVert="<<numOfVert<<", type="<<type);

  
  if (type==0){
    size_tt base=2;
    adjList=new size_tt*[numOfVert+1];
    adjList[0]=new size_tt[numOfVert];
    adjList[0][0]=base;
    adjList[1] = new size_tt[h-1];
    
    for(int i=1;i<first_leaf;i++) {
      adjList[0][i]=base+1;
      adjList[i+1] = new size_tt[h-1];
    }
    for(int i=first_leaf;i<numOfVert;i++) {
      adjList[0][i]=1;
      adjList[i+1] = new size_tt[h-1];
    }
    
    //setting adjacencies for the root  
    for (int i=0;i<base;i++)
      adjList[1][i]=i+1;
    
    //setting adjacencies for the internal nodes
    
    for(int i=1;i<first_leaf;i++){
      adjList[i+1][0]=(size_tt)((i-1)/base);
      for(int j=1;j<=base;j++)
	adjList[i+1][j]=base*i+j;
    }
    
    //setting adjacencies for the leafs
    for(int i=first_leaf;i<numOfVert;i++) 
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
	size_tt th=rand()%(numOfVert*numOfVert);
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
//      elements (we can set 'len' to be any number <= lenght[array])
//      put the first newLen elements of 'array' size_tto newArray.
//
//****************************************************************
 void Graph::rand_Perm(size_tt *array,
                      size_tt len,
                      size_tt *newArray,
                      size_tt newLen)
{
//    sfast_Rand(time(NULL));
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
//	Method name : (<<)
//
//	Description : print the compact adjacency list of the graph
//
//****************************************************************
ostream &operator<<( ostream &output,
                     const Graph &graph )
{
#if DEBUG
    debug("Entering (<<) " );
#endif

    for(size_tt index1 = 0; index1 < graph.numOfVert ; index1++)
    {
        output << index1 << ": ";
        
        for(size_tt index2 = 0;
            index2 < graph.adjList[0][index1]; index2++)
        {
            output << graph.adjList[ index1 + 1][ index2 ] << ' ';
        }
        output << endl;
    }
//    output );

#if DEBUG
    debug("Leaving (<<) " );
#endif

    return output;
}

//****************************************************************
//
//	Method name : print_Cpt_Adj_List_To_File
//
//****************************************************************
void Graph::print_Cpt_Adj_List_To_File() 
{
#if DEBUG
    debug("Entering print_Cpt_Adj_List_To_File" );
#endif
    
    int history;
    ifstream inHistFile("Graphs/history", ios::in);
    if(!inHistFile)
    {
        cerr << "Cannot open Graphs/history";
        exit(1);
    }
    inHistFile >> history;
    inHistFile.close();
    ofstream outHistFile("Graphs/history", ios::out);
    if(!outHistFile)
    {
        cerr << "Cannot write to Graphs/history";
        exit(1);
    }
    outHistFile << (history + 1);
    
    // create a name in the format Graph/numOfVert/GraphNumber
    string Graphs = "Graphs/";
    string num = itoa(numOfVert);
    string hyphen = "_";
    string historyStr = itoa(history);
    string fileName = Graphs + num + hyphen  + historyStr;
    const char * file = fileName.c_str();
    ofstream out(file);
    if(!out)
    {
        debug("File " << file << " could not be opened");
        exit(1);
    }
    else
        out << *this;

#if DEBUG
    debug("Exiting print_Adj_List_To_File" );
#endif
}

//****************************************************************
//
//	Method name : read_Cpt_Adj_List_From_File
//
//****************************************************************
void Graph::read_Cpt_Adj_List_From_File(const char *file) 
{
#if DEBUG
    debug("Entering read_Adj_List_From_File");
#endif

    ifstream in(file, ios::in);
    if(!in)
    {
        cerr << "Cannot open " << file;
        exit(1);
    }

    char c;
    size_tt index;
    
    in >> numOfVert;

    // initializing adjList
    adjList   = new size_tt*[ numOfVert ];
    adjList[0] = new size_tt[numOfVert - 1];

    for(size_tt adjVert = 0; adjVert < numOfVert - 1; adjVert++)
        in >> adjList[ 0 ][ adjVert ];

    
    for(size_tt vert = 1; vert < numOfVert; vert++) {
        in >> index; //read but do not record the "(vert #):"
        in >> c;
        adjList[ vert ] = new size_tt[adjList[ 0 ][ vert-1 ]];
        
        for(size_tt adjVert = 0; adjVert < adjList[ 0 ][ vert-1 ];
            adjVert++){
            in >> adjList[ vert ][ adjVert ];
        }
    }

#if DEBUG
    debug("Exiting read_Cpt_Adj_List_From_File");
#endif
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
//    return maximal degree vertex, its degree, and also
//    minimal degree of the graph
//
//****************************************************************
size_tt Graph::get_Max_Deg_Vert(size_tt &max_deg, size_tt &min_deg) const
{
    size_tt max_deg_vert = 0;
    max_deg = 0;
    min_deg = numOfVert;
    
#if DEBUG_GRAPH
    debug("numOfVert=" << numOfVert);
    debug("adjList:\n" << adjList );
#endif        
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
//	Method name : read_Graph_From_IWB_File
//
//      reading (Frick's) Tom Sawyer format
//      use temp array of size numOfVert
//      place the result in the compact adjList
//      and then translate it to full format
//
//**************************************************************
void Graph::read_Graph_From_IWB_File(const char *file)
{
#if DEBUG
    debug("ENTERING");
#endif

    size_tt _numOfEdges;
    size_tt v1, v2;
    char c;
    char str[10];

    ifstream in(file);
    if(!in){
        debug("Cannot open the file" << file );
        exit(1);
    }

    in >> numOfVert >> _numOfEdges >> c;
    
//      debug("numOfVert=" << numOfVert
//            << ", _numOfEdges=" << _numOfEdges
//            << ", c=" << c);
    
    // read a list of vertices, but don't record it
    for(size_tt index = 0; index < numOfVert; index++){
        in >> v1 >> str >> str;
//            debug("v1=" << v1 << ", str=" << str);
    }
    
    int p = in.tellg();
    
    size_tt *degs = new size_tt[numOfVert];
    
    for(size_tt index = 0; index < _numOfEdges; index++){
        in >> v1 >> v2;
//        debug("v1=" << v1 << ", v2=" << v2);
        degs[v1-1] += 1;
        degs[v2-1] += 1;            
    }
    in.seekg(p, ios::beg);
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];
    adjList[0] = degs;
    for(size_tt vert = 1; vert < numOfVert+1; vert++)
        adjList[ vert ] = new size_tt[degs[vert-1]];
    
    size_tt *count = new size_tt[numOfVert];
    for(size_tt i = 0; i < numOfVert; i++)
        count[i] = 0;
    
    for(size_tt index = 0; index < _numOfEdges; index++){
        in >> v1 >> v2;
//            debug("2nd time: v1=" << v1 << ", v2=" << v2);
        adjList[v1][count[v1-1]++] = v2-1;
        adjList[v2][count[v2-1]++] = v1-1;            
    }
//    debug("adjList: " << adjList);
    
        in.close();
	delete [] count;
        //delete [] degs;

#if DEBUG
    debug("LEAVING");
#endif
}

//**************************************************************
//
//	square_Torus()
//
//**************************************************************
void Graph::square_Torus( size_tt base )
{
#if DEBUG
    debug("Entering");
#endif 

    cycle_Graph( base );
    for(size_tt i = 0; i < 2; i++)
        double_Graph();
    
#if DEBUG
    debug("Leaving");
#endif 
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
  for(int i=0;i<depth;i++)
    total = total + (size_tt)pow(base,i);
  numOfVert = total;

  adjList=new size_tt*[numOfVert+1];
  adjList[0]=new size_tt[numOfVert];
  adjList[0][0]=base;
  adjList[1] = new size_tt[base];
  
  int first_leaf=numOfVert-(size_tt)pow(base,depth-1);

  for(int i=1;i<first_leaf;i++) {
    adjList[0][i]=base+1;
    adjList[i+1] = new size_tt[base+1];
  }
  for(int i=first_leaf;i<numOfVert;i++) {
    adjList[0][i]=1;
    adjList[i+1] = new size_tt[1];
  }

  //setting adjacencies for the root  
  for (int i=0;i<base;i++)
    adjList[1][i]=i+1;

  //setting adjacencies for the internal nodes

  for(int i=1;i<first_leaf;i++){
    adjList[i+1][0]=(size_tt)((i-1)/base);
    for(int j=1;j<=base;j++)
      adjList[i+1][j]=base*i+j;
  }
		 
  //setting adjacencies for the leafs
  for(int i=first_leaf;i<numOfVert;i++) 
    adjList[i+1][0]=(size_tt)((i-1)/base);


//    debug("adjList: ");
//    for(size_tt v=0; v< numOfVert; v++){
//      cout << v << ": (adjList[0]["<<v<<"]="<<adjList[0][v]<<") : ";
//      for(size_tt j=0; j < adjList[0][v]; j++){
//        cout << adjList[v+1][j] << ' ';
//      }
//      cout << endl;
//    }
//    exit(1);
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
    debug("h=" << h<<", w="<< w<<", numOfVert=" << numOfVert );
    
    adjList = new size_tt*[numOfVert+1];
    adjList[0] = new size_tt[numOfVert];

    //setting degrees of the first colunm
    adjList[0][0]=2;
    adjList[0][h-1]=2;
    for(int k=1; k < h-1; k++)
        adjList[0][k]=3;

//      for(int v=0; v < h; v++)
//          debug("adjList[0]["<<v<<"]="<<adjList[0][v]);
    
    // memory alloc for the adj list of the elements of the first column
    for(size_tt v = 0; v < h; v++)
        adjList[v+1] = new size_tt[adjList[0][v]];

    // setting degrees and memory alloc of the elements of the remaining columns
    for(int i=1; i<w; i++){     
//        debug("i=" << i);
        if( i == w-1 ){
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v];
                adjList[ (v + i*h) + 1 ] = new size_tt[adjList[0][v + i*h]];
//                debug("adjList[0]["<<v+i*h<<"]="<<adjList[0][v + i*h]);
            }
        } else {
            for(size_tt v = 0; v < h; v++){
                adjList[ 0 ][ v + i*h ] = adjList[0][v]+1;
                adjList[ v+1+ i*h ] = new size_tt[adjList[0][ v + i*h]];
//                debug("adjList[0]["<<v+i*h<<"]="<<adjList[0][ v + i*h ]); 
	    }
	}
    }
    //
    //  SETTING ADJACENCY LISTS OF ALL VICES
    //

    // interior vertices
    for(int i=1; i < w-1; i++)
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
    for(int i=1; i < w-1; i++){
        adjList[ i*h + 1][0] = (i-1)*h;
        adjList[ i*h + 1][1] = 1 + i*h;
        adjList[ i*h + 1][2] = (i+1)*h;
    }

    // bottom horizontal edge
    for(int i=1; i < w-1; i++){
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
    
//      debug("adjList: ");
//      for(size_tt v=0; v< numOfVert; v++){
//          cout << v << ": (adjList[0]["<<v<<"]="<<adjList[0][v]<<") : ";
//          for(size_tt j=0; j < adjList[0][v]; j++){
//              cout << adjList[v+1][j] << ' ';
//          }
//          cout << endl;
//      }
//      exit(1);
    

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


//       debug("adjList: ");
//       for(size_tt v=0; v< numOfVert; v++){
//           cout << v << ": (adjList[0]["<<v<<"]="<<adjList[0][v]<<") : ";
//           for(size_tt j=0; j < adjList[0][v]; j++){
//               cout << adjList[v+1][j] << ' ';
//          }
//           cout << endl;
//       }
      //      exit(1);
    

#if DEBUG_MISH_GRAPH 
    debug("Leaving mesh_Graph");
#endif
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

  for(int i=2;i<=h;i++) {
    for(int j=1;j<=i;j++){
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
//        debug("adjList: ");
//        for(size_tt v=0; v< numOfVert; v++){
//            cout << v << ": (adjList[0]["<<v<<"]="<<adjList[0][v]<<") : ";
//            for(size_tt j=0; j < adjList[0][v]; j++){
//                cout << adjList[v+1][j] << ' ';
//           }
//            cout << endl;
//        }
       //       exit(1);
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

  for(int i=2;i<=h;i++) {
    for(int j=1;j<=i;j++){
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
	  adjList[index+1]=new size_tt[2];
	  adjList[index+1][0]=index-(i-1);
	  adjList[index+1][1]=index+1;
	  adjList[index+1][2]=numOfVert-1;
	}
	else if(j==i) {
	  adjList[0][index]=3;
	  adjList[index+1]=new size_tt[2];
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
    for(int i=1;i<=h;i++)
      temp=temp+(size_tt)pow(3,i);
    numOfVert=temp;
    debug("numOfVert="<<numOfVert);
    adjList = new size_tt*[numOfVert+1];
    adjList[0]=new size_tt[numOfVert];
    for(int i=0;i<numOfVert;i++) {
      adjList[i+1] = new size_tt[4];
      adjList[0][i]=0;
    }
    
    sierpinski_recurse(h,0,0,1,2);
  }
  else {
    numOfVert = sierpinski_recurse3D(h,0,0,1,2,3,false);
//    debug("numOfVert="<<numOfVert);
    adjList = new size_tt*[numOfVert+1];
    adjList[0]=new size_tt[numOfVert];
    for(int i=0;i<numOfVert;i++) {
      adjList[i+1] = new size_tt[6];
      adjList[0][i]=0;
    }
    vert=4;
    sierpinski_recurse3D(h,0,0,1,2,3, true);
  }
}

//**************************************************************
//
//	itoa()
//
//	convert an int into a string object
//
//**************************************************************
string Graph::itoa(int n) 
{
    int m = n;
    int index = 0;
    int base = 10;
    int sign;
    
    // determine the number of digits in n
    int numOfDig = 1;
    while (( m /= 10) > 0)
        numOfDig++;

    m = n;
    char *cstrg = new char[numOfDig+1];

    if ((sign = m) < 0)     // record sign 
        m = -(m + 1);       // make n positive but then subtract 1, to
                            // make sure it does not go over the max_int
    do
    {                       // generate digits in reverse order 
        if (sign < 0 && index == 0)
            cstrg[index++] = m % base + '1';
        else
            cstrg[index++] = m % base + '0';   // get next digit
    }
    while (( m /= base) > 0);
    
    if( sign < 0)
        cstrg[index++] = '-';
    
    cstrg[index] = '\0';

    reverse(cstrg);
    string strg(cstrg);

    delete [] cstrg;
    return strg;
}

//**************************************************************
//
//	Method name : reverse
//
//	Description : reverse string s in place 
//
//**************************************************************
void Graph::reverse(char s[]) 
{
    int c, i, j;
    
    for (i = 0, j = strlen(s) - 1; i < j; i++, j--)
        c = s[i], s[i] = s[j], s[j] = c;
}

