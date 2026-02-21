//MesaPlot.cpp - member function definitions for MesaPlot class.

#include "MesaPlot.hpp"

#define DEBUG_MESA 0

#define ENG_MISH_v5          12
#define ENG_MISH_v6          13

bool mlistSwitch = true;

DrawGraph *MesaPlot::dgPtr = 0;
Graph const *MesaPlot::graph = 0;
int MesaPlot::dim = 2;
size_tt MesaPlot::slowTime = 0;
size_tt MesaPlot::engf = 6;
int MesaPlot::numOfVert = 1;
int MesaPlot::width = 1000; 
int MesaPlot::height = 1000;
GLuint MesaPlot::blueMarble = 0;
GLuint MesaPlot::blueMarble2 = 1;
GLuint MesaPlot::redMarble2 = 2;
GLuint MesaPlot::bluePoint = 3;
GLuint MesaPlot::redPoint = 4;
GLuint MesaPlot::yellowPoint = 5;
bool MesaPlot::createList = false;
int MesaPlot::xspin = 0;
int MesaPlot::yspin = 0;
int MesaPlot::zspin = 0;
const int MesaPlot::VIEW_TURN_RATE = 2;
const float MesaPlot::ZOOM_FACTOR = .1;
Point<> MesaPlot::baricenter;
Point4D<> MesaPlot::baricenter4D;
bool MesaPlot::newPos = false;
float MesaPlot::dist = 1;
GLuint MesaPlot::finalDrawing;
coord_t MesaPlot::boxSize = 1;
bool MesaPlot::color = false;
Point4D<GLfloat> MesaPlot::projVect;
Point4D<GLfloat> MesaPlot::nProjVect;
Point4D<GLfloat> MesaPlot::e1;
Point4D<GLfloat> MesaPlot::e2;
Point4D<GLfloat> MesaPlot::e3;

//**************************************************************
//
//      class constructor
//
//**************************************************************
MesaPlot::MesaPlot( const Graph *_graph,
                    DrawGraph *e       ,
                    size_tt _engf      ,
                    size_tt _slow      ,
                    int _width         ,
                    int _height        ,
                    bool _color)
{
#define DEBUGa 0
#if DEBUGa
    cout << "Entering MesaPlot constructor1" << endl;
#endif
    graph         = _graph;
    dgPtr         = e;
    dim           = e->get_Dim();
    engf          = _engf;
    slowTime      = _slow;
    numOfVert     = e->get_NumOfVert();
    color         = _color;
    width         = _width;
    height        = _height;
    boxSize       = (coord_t)(e->get_Edge() * .7 * e->get_Diam());
    //coord_t box2Size = 2 * boxSize + 1;
    int box2Size = 2 * boxSize + 1;

    if( dim == 4 ){
    srand(time(NULL));
        
        projVect  = Point4D<GLfloat>
            ((GLfloat)(rand() % box2Size) - boxSize,
             (GLfloat)(rand() % box2Size) - boxSize,
             (GLfloat)(rand() % box2Size) - boxSize,
             (GLfloat)(rand() % box2Size) - boxSize);

        nProjVect = projVect/projVect.fnorm();

        //
        // COMPUTING ORTHONORMAL BASIS OF THE PLANE PERPENDICULAR TO PROJVECT
        //
        // choosing e1
        e1 = Point4D<GLfloat>((GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize);
    
        float pe1 = nProjVect.getX() * e1.getX() +
            nProjVect.getY() * e1.getY() +
            nProjVect.getZ() * e1.getZ() +
            nProjVect.getW() * e1.getW();
        e1 -= pe1 * nProjVect;
        e1 /= (double)e1.fnorm();

        // choosing e2
        e2 = Point4D<GLfloat>((GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize);

        float pe2 = nProjVect * e2;
        float e1e2 = e1 * e2;
    
        e2 -= pe2 * nProjVect;
        e2 -= e1e2 * e1;
        e2 /= (double)e2.fnorm();

        while( (0.9 < e1e2 && e1e2 < 1.1) || (-1.1 < e1e2 && e1e2 < -0.9)){
            e2 = Point4D<GLfloat>((GLfloat)(rand() % box2Size) - boxSize,
                                  (GLfloat)(rand() % box2Size) - boxSize,
                                  (GLfloat)(rand() % box2Size) - boxSize,
                                  (GLfloat)(rand() % box2Size) - boxSize);
        
            pe2 = nProjVect.getX() * e2.getX() +
                nProjVect.getY() * e2.getY() +
                nProjVect.getZ() * e2.getZ() +
                nProjVect.getW() * e2.getW();
            e2 -= pe2 * nProjVect;
            e1e2 = e1 * e2;
            e2 -= e1e2 * e1;
            e2 /= (double)e2.fnorm();
        }
    
        // choosing e3
        e3 = Point4D<GLfloat>((GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize,
                              (GLfloat)(rand() % box2Size) - boxSize);

        float pe3 = nProjVect.getX() * e3.getX() +
            nProjVect.getY() * e3.getY() +
            nProjVect.getZ() * e3.getZ() +
            nProjVect.getW() * e3.getW();
        e3 -= pe3 * nProjVect;
        float e1e3 = e1 * e3;
        e3 -= e1e3 * e1;
        float e2e3 = e2 * e3;
        e3 -= e2e3 * e2;
        e3 /= (double)e3.fnorm();
        //checking if e3 is not too close to e1 or e2

        while( (0.9 < e1e3 && e1e3 < 1.1) || (-1.1 < e1e3 && e1e3 < -0.9) ||
               (0.9 < e2e3 && e2e3 < 1.1) || (-1.1 < e2e3 && e2e3 < -0.9)){
            e3 = Point4D<GLfloat>((GLfloat)(rand() % box2Size) - boxSize,
                                  (GLfloat)(rand() % box2Size) - boxSize,
                                  (GLfloat)(rand() % box2Size) - boxSize,
                                  (GLfloat)(rand() % box2Size) - boxSize);
        
            pe3 = nProjVect.getX() * e3.getX() +
                nProjVect.getY() * e3.getY() +
                nProjVect.getZ() * e3.getZ() +
                nProjVect.getW() * e3.getW();
            e3 -= pe3 * nProjVect;
            e1e3 = e1 * e3;
            e3 -= e1e3 * e1;
            e2e3 = e2 * e3;
            e3 -= e2e3 * e2;
            e3 /= (double)e3.fnorm();
        }

//      debug("projVect="<< projVect<<"\nnProjVect="<< nProjVect<<
//            "\ne1="<< e1<<", e1*nProjVect="<<e1*nProjVect<< 
//            "\ne2="<< e2<<", e2*nProjVect="<<e2*nProjVect<<", e1*e2="<<e1*e2<<
//            "\ne3="<< e3<<", e3*nProjVect="<<e3*nProjVect<<
//            ", e1*e3="<<e1*e3<<", e2*e3="<<e2*e3);
//      exit(1);
    }
#if DEBUGa
    cout << "Leaving MesaPlot constructor1" << endl;
#endif
}

//**************************************************************
//
//	Method name : lighting
//
//	Description : 
//
//**************************************************************
void MesaPlot::lighting()
{
    float x = 2 * width;
    
    GLfloat position[] = {x, x, x, 1.0};
    GLfloat ambient[] = { 1.0, 1.0, 1.0, 1.0 };
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glDepthFunc(GL_LESS);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    tbMatrix();
}

//**************************************************************
//
//	Method name : init
//
//	Description : 
//
//**************************************************************
void MesaPlot::init()
{
#if DEBUG_MESA
    cout << "Entering init" << endl;
#endif
   
   GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
   GLfloat mat_diffuse_blue[] = { 0.1, 0.5, 0.8, 1.0 };
   GLfloat mat_diffuse_red[] = { 1.0, 0.1, 0.1, 1.0 };
   GLfloat mat_amb_diff_yellow[]  = {1.0, 1.0, 0.0, 1.0};
   GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat high_shininess[] = { 100.0 };

   glMaterialfv(GL_FRONT, GL_AMBIENT, no_mat);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_blue);
   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
   glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);

   glEnable(GL_DEPTH_TEST);
   glShadeModel(GL_SMOOTH);

   tbInit(GLUT_LEFT_BUTTON);
  
   GLint shorterEdge = (width > height) ? height : width;
   GLfloat blueMarbleRadius;
   blueMarbleRadius = 10;

//   blueMarble2
   glNewList (blueMarble2, GL_COMPILE);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_amb_diff_yellow);
   glutSolidSphere(blueMarbleRadius * boxSize / shorterEdge, 5, 5);
   glEndList();

//   redMarble2 
   glNewList (redMarble2, GL_COMPILE);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse_red);
   glutSolidSphere(blueMarbleRadius * boxSize /shorterEdge, 5, 5);
   glEndList();
   
#if DEBUG_MESA
    cout << "Exiting init" << endl;
#endif
}

void MesaPlot::init2()
{
#if DEBUG_MESA
    cout << "Entering init2" << endl;
#endif

    GLfloat mat_amb_diff_red[]   = {1.0, 0.0, 0.0, 1.0};
    GLfloat mat_amb_diff_blue[]    = {0.0, 0.0, 1.0, 1.0};
    GLfloat mat_amb_diff_yellow[]  = {1.0, 1.0, 0.0, 1.0};

   if(color)
       glClearColor(1.0, 1.0, 1.0, 0.0);
   else
       glClearColor(0.0, 0.0, 0.0, 0.0);
   
//   bluePoint
   glNewList (bluePoint, GL_COMPILE);
   glPointSize(1.5);
   glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_diff_blue);
   glBegin(GL_POINTS);
   glVertex3f(0,0,0);
   glEnd();
   glEndList();

// yellowPoint
   glNewList (yellowPoint, GL_COMPILE);
   glPointSize(1.5);
   glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_diff_yellow);
   glBegin(GL_POINTS);
   glVertex3f(0,0,0);
   glEnd();
   glEndList();

//   redPoint
   glNewList (redPoint, GL_COMPILE);
   glPointSize(1.5);
   glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb_diff_red);
   glBegin(GL_POINTS);
   glVertex3f(0,0,0);
   glEnd();
   glEndList();

#if DEBUG_MESA
    cout << "Exiting init2" << endl;
#endif
}



//**************************************************************
//
//	Method name : reshape
//
//	Description : 
//
//**************************************************************
 void MesaPlot::reshape(int w, int h)
{
#if DEBUG_MESA
    cout << "Entering reshape" << endl;
#endif
    GLdouble margin = 0.1;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    tbReshape(w, h);
            
    if (w <= h ) {
        glOrtho (-(1.0 + margin) * boxSize,
                 (1.0 + margin) * boxSize,
                 -(((GLfloat)h)/(GLfloat)w + margin) * boxSize,
                 (((GLfloat)h)/(GLfloat)w + margin) * boxSize,
                 -10.0 * boxSize, 10.0 * boxSize);
    }
    else
        glOrtho (-margin,(GLfloat)w/((GLfloat)h) + margin,
                 -margin, 1.0 + margin,
                 -10.0, 10.0);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

#if DEBUG_MESA
    cout << "Exiting reshape" << endl;
#endif
}

//**************************************************************
//
//	Method name : keyboard
//
//	Description : 
//
//**************************************************************
 void MesaPlot::keyboard(unsigned char key, int x, int y)
{
    switch (key) {
        case 'r':
            glutIdleFunc(MesaPlot::engine);
            break;
        case 's':
            glutIdleFunc(NULL);
            break;
        case 'f':
            move_Forward();
            break;
        case 'b':
            move_Backward();
            break;
        case 'x':
            glutIdleFunc(MesaPlot::x_Rotate_CC);
            break;
        case 'X':
            glutIdleFunc(MesaPlot::x_Rotate_CW);
            break;
        case 'y':
            glutIdleFunc(MesaPlot::y_Rotate_CC);
            break;
        case 'Y':
            glutIdleFunc(MesaPlot::y_Rotate_CW);
            break;
        case 'z':
            glutIdleFunc(MesaPlot::z_Rotate_CC);
            break;
        case 'Z':
            glutIdleFunc(MesaPlot::z_Rotate_CW);
            break;
        case 27:
        case 'q':
        case 'Q':
            exit(0);
            break;
            
    }
}

void MesaPlot::special(int key, int x, int y)
{
    switch (key){
        case GLUT_KEY_RIGHT:
                glutIdleFunc(MesaPlot::z_Rotate_CW);
            break;
        case GLUT_KEY_LEFT:
                glutIdleFunc(MesaPlot::z_Rotate_CC);
            break;
        case GLUT_KEY_UP:
                glutIdleFunc(MesaPlot::x_Rotate_CW);
            break;
        case GLUT_KEY_DOWN:
                glutIdleFunc(MesaPlot::x_Rotate_CC);
            break;
    }
}

//**************************************************************
//
//	Method name : mouse
//
//	Description : 
//
//**************************************************************
void MesaPlot::mouse(int button, int state, int x, int y) 
{
   switch (button) { 
      case GLUT_MIDDLE_BUTTON:
         if (state == GLUT_DOWN)
            glutIdleFunc(MesaPlot::engine);
         break;
      case GLUT_LEFT_BUTTON:
          tbMouse(button, state, x, y);
          break;
      case GLUT_RIGHT_BUTTON:
         if (state == GLUT_DOWN)
            glutIdleFunc(NULL);
         break;
      default:
         break;
   }
}
   
void MesaPlot::motion(int x, int y){
  tbMotion(x, y);
}

//**************************************************************
//
//	Method name : draw_Graph
//
//	Description : 
//
//**************************************************************
void MesaPlot::draw_Graph(int argc,
                          char** argv)
{
#if DEBUG_MESA
    cout << "Entering draw_Graph" << endl;
#endif

    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize (width, height);
    glutInitWindowPosition(205, 25);
    glutCreateWindow(argv[0]);
    init();
    init2();
    glutReshapeFunc(MesaPlot::reshape);
    glutDisplayFunc(MesaPlot::display);
    glutKeyboardFunc (MesaPlot::keyboard);
    glutSpecialFunc(MesaPlot::special);
    glutMotionFunc(MesaPlot::motion);
    glutMouseFunc(MesaPlot::mouse);
    glutIdleFunc(MesaPlot::engine);

    glutCreateMenu(menu);
    glutAddMenuEntry("Write picture1.tiff", 1);
    glutAddMenuEntry("Write picture2.tiff", 2);
    glutAddMenuEntry("Write picture3.tiff", 3);
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutMainLoop();

#if DEBUG_MESA
    cout << "Exiting draw_Graph" << endl;
#endif
}

//**************************************************************
//
//	Method name : engine
//
//	Description : 
//
//**************************************************************
void MesaPlot::engine()
{
    static bool localSwitch = true;
    static bool turnIdleOff = false;

//      Graph g;
//      string time = g.itoa(slowTime);
//      string sleep = "sleep ";
//      string sleepT = sleep + time;
//      const char * command = sleepT.c_str();
        
      if(slowTime)
          system("sleep 1s");

    createList = dgPtr->createList;
    
    if(engf == ENG_MISH_v5)
        dgPtr->mish_engine_v5();
    else if(engf == ENG_MISH_v6)
        dgPtr->mish_engine_v6();
    else if( localSwitch ) {
        createList = true;
        localSwitch = false;
    } else 
        createList = false;

    if(turnIdleOff)
        glutIdleFunc(NULL);
        
    if(createList)
        turnIdleOff = true;
    
    glutPostRedisplay();
}

//**************************************************************
//
//	Method name : move_Forward/Backward
//
//	Description : routines used to zoom the light and viewer
//      pos
//
//**************************************************************
void MesaPlot::move_Forward()
{
    dist += ZOOM_FACTOR;
    glutPostRedisplay();
}

void MesaPlot::move_Backward()
{
    dist -= ZOOM_FACTOR;
    glutPostRedisplay();
}

//**************************************************************
//
//	Method name : ?_Rotate_CW(CC)
//
//	Description : rotation about the ?-axes in the clockwise (CW)
//      or counter-clockwise (CC) direction
//
//**************************************************************

void MesaPlot::x_Rotate_CC()
{
    xspin = (xspin + VIEW_TURN_RATE) % 360;
    glutPostRedisplay();
}

void MesaPlot::x_Rotate_CW()
{
    xspin = (xspin - VIEW_TURN_RATE) % 360;
    glutPostRedisplay();
}

void MesaPlot::y_Rotate_CC()
{
    yspin = (yspin + VIEW_TURN_RATE) % 360;
    glutPostRedisplay();
}

void MesaPlot::y_Rotate_CW()
{
    yspin = (yspin - VIEW_TURN_RATE) % 360;
    glutPostRedisplay();
}

void MesaPlot::z_Rotate_CC()
{
    zspin = (zspin + VIEW_TURN_RATE) % 360;
    glutPostRedisplay();
}

void MesaPlot::z_Rotate_CW()
{
    zspin = (zspin - VIEW_TURN_RATE) % 360;
    glutPostRedisplay();
}


//**************************************************************
//
//	Method name : writetiff()
//
//	Description : this routine has been lifted up from writetiff.c
//      code by Mark J. Kilgard
//      it saves a screen shot in a TIFF file picture[1-3].tif
//      writetiff uses Sam Leffler's libtiff library to write to
//      a TIFF an image grabbed from the screen with glReadPixels. 
//
//**************************************************************
int MesaPlot::writetiff(const char *filename, const char *description,
                        int x, int y, int width, int height,
                        int compression)
{
  TIFF *file;  
  GLubyte *image, *p;
  int i;

  file = TIFFOpen(filename, "w");
  if (file == NULL) {
    return 1;
  }
  image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

  /* OpenGL's default 4 byte pack alignment would leave extra bytes at the
     end of each image row so that each full row contained a number of bytes
     divisible by 4.  Ie, an RGB row with 3 pixels and 8-bit componets would
     be laid out like "RGBRGBRGBxxx" where the last three "xxx" bytes exist
     just to pad the row out to 12 bytes (12 is divisible by 4). To make sure
     the rows are packed as tight as possible (no row padding), set the pack
     alignment to 1. */
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
  TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
  TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
  TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(file, TIFFTAG_COMPRESSION, compression);
  TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, description);
  p = image;
  for (i = height - 1; i >= 0; i--) {
    if (TIFFWriteScanline(file, p, i, 0) < 0) {
      free(image);
      TIFFClose(file);
      return 1;
    }
    p += width * sizeof(GLubyte) * 3;
  }
  TIFFClose(file);
  return 0;
}

//**************************************************************
//
//	Method name : menu()
//
//**************************************************************
void MesaPlot::menu(int value)
{
    switch (value) {
        case 1:
            writetiff("picture1.tif", "OpenGL-rendered gears", 0, 0, width, height,
                      COMPRESSION_NONE);
            break;
        case 2:
            writetiff("picture2.tif", "OpenGL-rendered gears", 0, 0, width, height,
                      COMPRESSION_NONE);
            break;
        case 3:
            writetiff("picture3.tif", "OpenGL-rendered gears", 0, 0, width, height,
                      COMPRESSION_NONE);
            break;
  }
}

/**
 * @brief Updates the barycenter of the graph and adjusts vertex positions.
 *
 * This function calculates the barycenter of the current set of vertices in the graph
 * and shifts all vertex positions so that the barycenter is at the origin (0,0,0).
 * The calculation and adjustment are performed only if the number of vertices has changed
 * since the last update.
 *
 * @details The function handles both 3D and 4D coordinate systems, determined by the 'dim' member variable.
 * For 4D coordinates, it uses 'baricenter4D' and 'pos4D', while for 3D it uses 'baricenter' and 'pos'.
 *
 * The function uses the following member variables:
 * - dgPtr: Pointer to the graph data structure
 * - dim: Dimension of the coordinate system (3 or 4)
 * - baricenter: 3D barycenter (used when dim == 3)
 * - baricenter4D: 4D barycenter (used when dim == 4)
 *
 * @note This function modifies the positions of vertices in the graph (dgPtr->pos or dgPtr->pos4D).
 * @note The function uses a static variable 'mSize' to keep track of the previous number of vertices.
 *
 * @see draw_vertices()
 * @see draw_edges()
 */
void MesaPlot::update_barycenter()
{
    size_tt misfSize = dgPtr->misfSize[dgPtr->misfLevel];
    static size_tt mSize = dgPtr->numOfInitVert;

    if (misfSize == dgPtr->numOfInitVert || misfSize != mSize) {
        if (dim == 4) {
            baricenter4D.set_to_zero();
            for (int v = 0; v < misfSize; v++)
                baricenter4D += dgPtr->pos4D[dgPtr->mish[v]];
            baricenter4D /= (coord_t)misfSize;

            // Shift position vectors so that the barycenter is at (0,0,0)
            for (int v = 0; v < misfSize; v++)
                dgPtr->pos4D[dgPtr->mish[v]] -= baricenter4D;
        } else {
            baricenter.set_to_zero();
            for (int v = 0; v < misfSize; v++)
                baricenter += dgPtr->pos[dgPtr->mish[v]];
            baricenter /= (coord_t)misfSize;

            // Shift position vectors so that the barycenter is at (0,0,0)
            for (int v = 0; v < misfSize; v++)
                dgPtr->pos[dgPtr->mish[v]] -= baricenter;
        }
        mSize = misfSize;
    }
}
/**
 * @brief Draws all vertices of the graph.
 *
 * This function iterates through all vertices in the current graph and renders them
 * using OpenGL commands. The appearance of each vertex (color and size) is determined
 * by its properties and the total number of vertices.
 *
 * @details The function handles both 3D and 4D coordinate systems, determined by the 'dim' member variable.
 * For 4D coordinates, it projects the vertex to 3D space using the 'project()' function before rendering.
 *
 * The function uses the following member variables:
 * - dgPtr: Pointer to the graph data structure
 * - dim: Dimension of the coordinate system (3 or 4)
 * - blueMarble2, yellowPoint, redMarble2, redPoint: Display lists for different vertex appearances
 *
 * @note This function assumes that the OpenGL state (including the model-view matrix) has been
 *       properly set up before it is called.
 * @note The function uses glPushMatrix() and glPopMatrix() for each vertex to isolate transformations.
 *
 * @see update_barycenter()
 * @see draw_edges()
 * @see project()
 */
void MesaPlot::draw_vertices()
{
    size_tt marbleTreshold = 80;
    size_tt misfSize = dgPtr->misfSize[dgPtr->misfLevel];

    for (int v = 0; v < misfSize; v++) {
        size_tt vert = dgPtr->mish[v];
        glPushMatrix();
        if (dim == 4) {
            GLfloat x, y, z;
            project(vert, x, y, z);
            glTranslatef(x, y, z);
        } else {
            glTranslatef((GLfloat)(dgPtr->pos[vert].getX()),
                         (GLfloat)(dgPtr->pos[vert].getY()),
                         (GLfloat)(dgPtr->pos[vert].getZ()));
        }

        if (dgPtr->inv[vert] >= dgPtr->prevMishSize) {
            if (misfSize < marbleTreshold)
                glCallList(blueMarble2);
            else
                glCallList(yellowPoint);
        } else {
            if (misfSize < marbleTreshold)
                glCallList(redMarble2);
            else
                glCallList(redPoint);
        }
        glPopMatrix();
    }
}

/**
 * @brief Draws all edges of the graph.
 *
 * This function iterates through all edges in the current graph and renders them
 * using OpenGL commands. The appearance of each edge (color) is determined by the
 * 'color' member variable.
 *
 * @details The function handles both 3D and 4D coordinate systems, determined by the 'dim' member variable.
 * For 4D coordinates, it projects the vertex endpoints to 3D space using the 'project()' function before rendering.
 *
 * The function uses the following member variables:
 * - dgPtr: Pointer to the graph data structure
 * - graph: Pointer to the graph adjacency list
 * - dim: Dimension of the coordinate system (3 or 4)
 * - color: Boolean determining the color of the edges
 * - numOfVert: Total number of vertices in the graph
 *
 * @note This function only draws edges if dgPtr->misfLevel is 0.
 * @note This function assumes that the OpenGL state (including the model-view matrix) has been
 *       properly set up before it is called.
 * @note The function uses glPushMatrix() and glPopMatrix() for each edge to isolate transformations.
 *
 * @warning This function does not perform bounds checking on array accesses. Ensure that
 *          numOfVert and dgPtr->deg[vert] are correctly set to avoid out-of-bounds access.
 *
 * @see update_barycenter()
 * @see draw_vertices()
 * @see project()
 */
void MesaPlot::draw_edges()
{
    if (dgPtr->misfLevel != 0) return;

    GLfloat mat_amb_diff_black[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat mat_amb_diff_green[] = {0.0, 1.0, 0.0, 1.0};

    for (int vert = 0; vert < numOfVert-1; vert++) {
        size_tt deg = dgPtr->deg[vert];
        for (size_tt adjVert = 0; adjVert < deg; adjVert++) {
            int av = graph->adjList[vert+1][adjVert];
            if (av > vert) {
                if (color)
                    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_amb_diff_black);
                else
                    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_amb_diff_green);
                glPushMatrix();
                glBegin(GL_LINES);
                if (dim == 4) {
                    GLfloat x, y, z;
                    project(vert, x, y, z);
                    glVertex3f(x, y, z);
                    project(av, x, y, z);
                    glVertex3f(x, y, z);
                } else {
                    glVertex3f((GLfloat)(dgPtr->pos[vert].getX()),
                               (GLfloat)(dgPtr->pos[vert].getY()),
                               (GLfloat)(dgPtr->pos[vert].getZ()));
                    glVertex3f((GLfloat)(dgPtr->pos[av].getX()),
                               (GLfloat)(dgPtr->pos[av].getY()),
                               (GLfloat)(dgPtr->pos[av].getZ()));
                }
                glEnd();
                glPopMatrix();
            }
        }
    }
}

/**
 * @function orig_display
 * @brief Renders a sequence of graphical approximations, culminating in a final static display using OpenGL.
 *
 * This function is responsible for rendering a series of graph approximations onto the screen using OpenGL. It
 * handles the creation of display lists, which are used to efficiently render the graphical content across
 * multiple frames. The function ensures that the final drawing remains displayed on the screen, without being
 * erased or overwritten by subsequent redraw operations.
 *
 * Detailed Description:
 * - The function begins by clearing the screen buffer and enabling lighting effects.
 * - If `createList` is true, a new display list (`finalDrawing`) is generated and compiled. This list will
 *   store the graphical content, including vertices and edges, with anti-aliasing enabled for smoother visuals.
 * - If `mlistSwitch` is true, the display list is compiled with the updated graph data (vertices and edges).
 *   The `update_barycenter` function is called to adjust vertex positions before drawing the graph.
 * - Once the list compilation is complete, `mlistSwitch` and `createList` are set to false to prevent further
 *   recompilation and unnecessary redraws.
 * - If `mlistSwitch` is false, the function uses the existing display list (`finalDrawing`) to render the
 *   final graph without modifying it. The graph is displayed with the applied transformations (rotation and scaling).
 * - The function concludes by swapping the display buffers to present the final rendered frame.
 *   To avoid clearing the final frame, `glutPostRedisplay()` is only called if the display list is still being
 *   created or updated.
 *
 * @note This function is designed to ensure that the final graphical output remains on the screen after all
 *       rendering operations are complete, addressing an issue where previous versions might have cleared the
 *       final drawing from the screen.
 *
 * @variables:
 * - `xspin`, `yspin`, `zspin`: Floating-point variables representing the rotation angles for the x, y, and z axes.
 * - `dist`: A floating-point variable controlling the scaling of the rendered graph.
 * - `createList`: A boolean flag indicating whether a new display list needs to be created.
 * - `mlistSwitch`: A boolean flag indicating whether the display list should be updated.
 * - `finalDrawing`: An integer holding the ID of the OpenGL display list used for rendering the graph.
 *
 * @dependencies:
 * - Requires OpenGL libraries for rendering (e.g., `gl`, `glut`, etc.).
 * - Assumes that the `lighting()`, `update_barycenter()`, `draw_vertices()`, and `draw_edges()` functions are
 *   implemented and properly handle the corresponding OpenGL operations.
 *
 * @usage:
 * - This function is typically used as the display callback in an OpenGL rendering loop, invoked by `glutMainLoop()`.
 * - Ensure that the necessary OpenGL context is initialized before calling this function.
 */
void MesaPlot::display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    lighting();

    if(createList){
        // draw nice (antialiased) edges
        finalDrawing = glGenLists(1);
        glEnable(GL_LINE_SMOOTH);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
    }

    if( mlistSwitch == true ){
        if(createList)
            glNewList(finalDrawing, GL_COMPILE);

        update_barycenter();
        draw_vertices();
        draw_edges();

        if( createList ) {
            glEndList();
            mlistSwitch = false;
            createList = false; // Disable further list creation
        }
    }
    else
    {   // draw finalDrawing
        lighting();
        glPushMatrix();
        glRotatef(xspin, 1.0, 0.0, 0.0);
        glRotatef(yspin, 0.0, 1.0, 0.0);
        glRotatef(zspin, 0.0, 0.0, 1.0);
        glScalef(dist, dist, dist);
        glCallList(finalDrawing);
        glPopMatrix();
    }

    glutSwapBuffers();

    // Only request a redraw if still creating the list
    if (createList || mlistSwitch) {
        glutPostRedisplay();
    }
}
