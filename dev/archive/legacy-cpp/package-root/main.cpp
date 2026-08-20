//      main.cpp file

//**************************************************************
//
//      a driver for DrawGraph package
//
//**************************************************************


#include <time.h>

#include "DrawGraph.hpp"
#include "MesaPlot.hpp"

#define DEBUG_MAIN 0

#define ENG_MISH_v5          12
#define ENG_MISH_v6          13

//global parameters
size_tt _rounds       = 20;
size_tt _finalRounds  = 10;
size_tt _tinit_factor = 6;
size_tt _p1           = 0;
size_tt _p2           = 0;
size_tt _numOfVert    = 5;
size_tt _numOfItr     = 1;
size_tt _dim          = 3;
size_tt _engf         = ENG_MISH_v5;
size_tt _slow         = 0;
float _roundsPar      = 1.0;
coord_t _edge         = 32;
bool _readFromFile    = false;
bool _printToFile     = false;
char fileName[30];
char _posFile[30];
std::string graphStr       = "complete";
size_tt _par          = 3;        // parameter for hypercube etc.
size_tt _displayPar   = 1;          // dispaly switch, 1 = ON
int _width            = 700;
int _height           = _width;
size_tt _numOfInitVert= 15;
size_tt _numOfNbrs    = 10;
float _r              = 0.15;  //parameters of update_Local_Temp_v3()
float _s              = 3.0;
size_tt _T            = 4;
bool _color           = false;

//routine processing prompt parameters
void cases(int argc, char** argv);

int main(int argc, char** argv)
{
   clock_t c0, c1;
   time_t t0, t1;

   t0 = time(NULL);
   c0 = clock();

   std::set_new_handler(DrawGraph::noMoreMemory);
    
   cases(argc, argv);
   Graph graph;

   if( graphStr == "complete" ){
       graph.complete_Graph( _numOfVert );
   } else if( graphStr == "random" )
       graph.rand_Graph( _numOfVert, _T );
   else if( graphStr == "path" )
       graph.path_Graph( _numOfVert );
   else if( graphStr == "cycle" )
       graph.cycle_Graph( _numOfVert );
   else if( graphStr == "hypercube" )
       graph.hyper_Cube( _numOfVert ); //one should test if
   // _numOfVert is a power of 2
   else if( graphStr == "cylinder" )
       graph.square_Cylinder( _numOfVert, _T);
   else if( graphStr == "torus" )
       graph.torus( _numOfVert, _T);
   else if( graphStr == "mesh" )
       graph.mesh( _numOfVert );
   else if (graphStr == "twistedtorus")
       graph.twistedTorus(_numOfVert, _T, _p1, _p2);
   else if (graphStr == "tree")
       graph.tree(_numOfVert, _T);
   else if (graphStr == "meshX")
       graph.meshX(_numOfVert);
   else if (graphStr == "meshT")
       graph.meshT(_numOfVert);
   else if (graphStr == "meshTX")
       graph.meshTX(_numOfVert);
   else if (graphStr == "sierpinski")
       graph.sierpinski(_numOfVert, _T);
   else
       graph.read_Graph_From_IWB_File(graphStr.c_str());

   DrawGraph dg( graph,
                 _dim,
                 _rounds,
                 _finalRounds,
                 _tinit_factor,
                 _engf,
                 _numOfInitVert,
                 _numOfNbrs ,
                 _r,  //parameters of update_Local_Temp_v3()
                 _s,
                 _displayPar);

//      if( _readFromFile )
//          dg.read_Positions_From_File(fileName);
    
   if( _displayPar == 1 ){
       MesaPlot mp(&graph, &dg, _engf, _slow, _width, _height, _color);
       mp.draw_Graph(argc, argv);
   } else if( _displayPar == 0 ){
       if(_engf == ENG_MISH_v5)
           dg.mish_engine_v5();
       else if(_engf == ENG_MISH_v6)
           dg.mish_engine_v6();
        
       c1 = clock();
       t1 = time(NULL);
        
       debug("\nUser time=" << (c1-c0)/(double)CLOCKS_PER_SEC
             <<"s\nReal time="<<difftime(t1,t0)<<"s");

//        exit(1);
        
       MesaPlot mp(&graph, &dg, 0, _slow, _width, _height, _color);
       mp.draw_Graph(argc, argv);
   }
    
   return 0;
}

void cases(int argc, char** argv){
    int c;
    for(int arg = 1; arg < argc; arg++){
        if( (argv[arg])[0] == '-' ){
            c = argv[arg][1];
            switch(c){
            case '#':
                _numOfVert = atoi(argv[++arg]);
                break;
            case 'r':
                _rounds = atoi(argv[++arg]);
                break;
            case 'R':
                _finalRounds = atoi(argv[++arg]);
                break;
	        case 'n':
	            _numOfItr = atoi(argv[++arg]);
                break;
            case 'i':
	            _numOfInitVert = atoi(argv[++arg]);
                break;
            case 'b':
	            _numOfNbrs = atoi(argv[++arg]);
                break;
            case 'c':
	            _color = (bool)atoi(argv[++arg]);
                break;
            case 'S':
	            _s = atof(argv[++arg]);
	            break;
            case 'T':
	            _T = atoi(argv[++arg]); //thickness for cylinder and torus
	            break;
            case 'd':
                _dim = atoi(argv[++arg]);
//                    _depth_ratio = atof(argv[++arg]);
                break;
            case 'k':
                if(!strcmp("mish_v5", argv[arg+1]))
                    _engf = ENG_MISH_v5;
                else if(!strcmp("mish_v6", argv[arg+1]))
                    _engf = ENG_MISH_v6;
                else
                    _engf = atoi(argv[++arg]);
                break;
            case 'e':
                _edge = atoi(argv[++arg]);
                break;
            case 'w':
                _width = atoi(argv[++arg]);
                break;
            case 't':
                _tinit_factor = atoi(argv[++arg]);
                break;
            case 'g': //type of graph: complete, path, cycle,
                //hypercube - needs extra parameter -P
                //to specify dim,
                //squareTorus - cycle * square - needs extra
                //parameter -P to specify the num of el in
                //cycle,
                //squareCylinder - path * square - needs extra
                //parameter -P to specify the num of el in
                //path,
                graphStr = argv[++arg];
                break;
            case 's':
//                    if( !strcmp(argv[++arg], "slow") )
                _slow = atoi(argv[++arg]);
                break;
            case 'P':  // extra parameter for hypercube etc.
                _par = atoi(argv[++arg]);
                break;
            case 'p':
                _printToFile = true;
                if( (argv[++arg])[0] != '-' )
                    strcpy(_posFile, argv[arg]);
                break;
            case 'D'://display ON/OFF switch, 0 = OFF, 1 = ON
                _displayPar = atoi(argv[++arg]);
                break;
            default:
                debug("illegal option\n" << (char)c);
                std::cout << "usage: main [#doiecsrgp]\n"
                          << "# numOfVert\n"
                          << "k kernel type\n"
                          << "r rounds\n"
                          << "t tinit_factor\n"
                          << "n numOfItr\n"
                          << "d depth_ratio\n"
                          << "e edge\n"
                          << "w width\n"
                          << "g graph type: complete, path, cycle, torus, cylinder, hypercube\n"
                          << "T extra parameter for torus, cylinder, and hypercube\n"
                          << std::endl;
            }
        }
    }
}

