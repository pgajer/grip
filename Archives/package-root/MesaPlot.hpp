// MesaPlot.hpp a header file for MesaPlot class

#ifndef MESA_PLOT
#define MESA_PLOT

#include <GLUT/glut.h>
#include <tiffio.h>
#include "DrawGraph.hpp"
#include "trackball.h"

typedef unsigned short size_tt;

//**************************************************************
//
//	Class name : MesaPlot
//
//**************************************************************
class MesaPlot
{
public:

    // class constructor
    MesaPlot( const Graph *_graph,
              DrawGraph *e      ,
              size_tt _engf     ,
              size_tt _slow  = 0,
              int _width  = 700 ,
              int _height = 700,
              bool _color = true);

    ~MesaPlot(){
        delete graph;
        delete dgPtr;
    }
    // we need a copy constructor and an assignment operator !!!
    
    static void draw_Graph(int argc, char** argv);
    
private:
    static const Graph *graph;
    static DrawGraph *dgPtr;
    static int dim;           // do we need it (take it from DrawGraph)
    static size_tt engf;      // engine flag
    static size_tt slowTime;  // display one frame every slowTime second(s)
    static int numOfVert;     // do we need it (take it from DrawGraph or Graph)
    static int width, height;  // width and height of the display window
    static GLuint blueMarble;
    static GLuint blueMarble2;
    static GLuint redMarble2;
    static GLuint bluePoint;
    static GLuint redPoint;
    static GLuint yellowPoint;
    static GLuint finalDrawing;
    static int xspin;
    static int yspin;
    static int zspin;
    static bool createList;
    static coord_t boxSize;
    static bool color;

    static Point4D<GLfloat> projVect;   // projection vector
    static Point4D<GLfloat> nProjVect;  // projVect/|projVect|
    static Point4D<GLfloat> e1;    // basis vectors of an orthonomla base
    static Point4D<GLfloat> e2;    // of the plane perpendicular (in R^4)
    static Point4D<GLfloat> e3;    // to projVect

    static bool newPos;
    static const int VIEW_TURN_RATE;
    static const float ZOOM_FACTOR;
    static Point<> baricenter; 
    static Point4D<> baricenter4D; 
    static float dist;

    // OpenGL parameters and routines
    static Point<> get_Baricenter();
    static Point4D<> get_Baricenter4D();
    static void engine();

    static void draw_vertices();
    static void draw_edges();
    static void update_barycenter();
    static void display();

    static void reshape(int width, int height);
    static void keyboard(unsigned char key, int x, int y);
    static void special(int key, int x, int y);
    static void mouse(int button, int state, int x, int y);
    static void init();
    static void init2();
    static void lighting();
    static void x_Rotate_CC();
    static void x_Rotate_CW();
    static void y_Rotate_CC();
    static void y_Rotate_CW();
    static void z_Rotate_CC();
    static void z_Rotate_CW();
    static void move_Forward();
    static void move_Backward();

    static int writetiff(const char *filename, const char *description,
                         int x, int y, int width, int height,
                         int compression);
    static void menu(int value);
    static void motion(int x, int y);
    static void objects();
    static void project(size_tt vert,
                        GLfloat &x, GLfloat &y, GLfloat &z)
        {
            // first we compute the projection of pos[vert] onto the plane p
            // perpendicular to projVect

            // scalar product of nProjVect and dgPtr->pos[vert]
            float pv =
                nProjVect.getX() * dgPtr->pos4D[vert].getX() +
                nProjVect.getY() * dgPtr->pos4D[vert].getY() +
                nProjVect.getZ() * dgPtr->pos4D[vert].getZ() +
                nProjVect.getW() * dgPtr->pos4D[vert].getW();

            Point4D<GLfloat> pvProjVect = pv * nProjVect;
            
            float proj_x = dgPtr->pos4D[vert].getX() - pvProjVect.getX();
            float proj_y = dgPtr->pos4D[vert].getY() - pvProjVect.getY();
            float proj_z = dgPtr->pos4D[vert].getZ() - pvProjVect.getZ();
            float proj_w = dgPtr->pos4D[vert].getW() - pvProjVect.getW();
//          Point4D<> proj =
//           dgPtr->pos4D[vert] - (nProjVect * dgPtr->pos4D[vert]) * nProjVect;

            // next we compute the coordinates of proj in an orthonormal base
            // e1, e2, e3 of the plane p (calculated in the constructor of
            // MesaPlot.cpp)
            x = proj_x * e1.getX() +
                proj_y * e1.getY() +
                proj_z * e1.getZ() +
                proj_w * e1.getW();
            y = proj_x * e2.getX() +
                proj_y * e2.getY() +
                proj_z * e2.getZ() +
                proj_w * e2.getW();
            z = proj_x * e3.getX() +
                proj_y * e3.getY() +
                proj_z * e3.getZ() +
                proj_w * e3.getW();
        }
};

#endif
