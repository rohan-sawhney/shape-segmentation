#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "Mesh.h"

int gridX = 600;
int gridY = 600;
int gridZ = 600;

const double fovy = 50.;
const double clipNear = .01;
const double clipFar = 1000.;
double x = 0;
double y = 0;
double z = -2.5;

std::string path = "/Users/rohansawhney/Desktop/developer/C++/shape-segmentation/bunny1.obj";
Mesh mesh;
bool success = true;
int clusters = 4;
std::vector<Eigen::Vector3d> colors;

void printInstructions()
{
    std::cerr << "' ': segment\n"
              << "→/←: increase/decrease clusters\n"
              << "↑/↓: move in/out\n"
              << "w/s: move up/down\n"
              << "a/d: move left/right\n"
              << "r: reload mesh\n"
              << "escape: exit program\n"
              << std::endl;
}

void generateClusterColors()
{
    colors.clear();
    for (int i = 0; i < 360; i += 360 / clusters) {
        Eigen::Vector3d c;
        c.x() = (double)rand() / RAND_MAX;
        c.y() = (double)rand() / RAND_MAX;
        c.z() = (double)rand() / RAND_MAX;
        
        colors.push_back(c);
    }
}

void init()
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable(GL_DEPTH_TEST);
}

void draw()
{
    glLineWidth(1.0);
    glBegin(GL_LINES);
    
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        
        if (f->isBoundary()) continue;
        Eigen::Vector3d &color(colors[f->cluster]);
        glColor4f(color.x(), color.y(), color.z(), 0.5);
        
        glBegin(GL_LINE_LOOP);
        HalfEdgeCIter he = f->he;
        do {
            glVertex3d(he->vertex->position.x(), he->vertex->position.y(), he->vertex->position.z());
            
            he = he->next;
            
        } while (he != f->he);
        
        glEnd();
    }
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double aspect = (double)viewport[2] / (double)viewport[3];
    gluPerspective(fovy, aspect, clipNear, clipFar);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    gluLookAt(0, 0, z, x, y, 0, 0, 1, 0);
    
    if (success) {
        draw();
    }

    glutSwapBuffers();
}

void keyboard(unsigned char key, int x0, int y0)
{
    switch (key) {
        case 27 :
            exit(0);
        case 'a':
            x -= 0.03;
            break;
        case 'd':
            x += 0.03;
            break;
        case 'w':
            y += 0.03;
            break;
        case 's':
            y -= 0.03;
            break;
        case 'r':
            success = mesh.read(path);
            break;
        case ' ':
            if (success) {
                generateClusterColors();
                mesh.segment(clusters);
            }
            break;
    }
    
    glutPostRedisplay();
}

void special(int i, int x0, int y0)
{
    switch (i) {
        case GLUT_KEY_UP:
            z += 0.03;
            break;
        case GLUT_KEY_DOWN:
            z -= 0.03;
            break;
        case GLUT_KEY_LEFT:
            if (clusters > 1) clusters--;
            break;
        case GLUT_KEY_RIGHT:
            clusters++;
            break;
    }
    
    std::string title = "Shape Segmentation, clusters: " + std::to_string(clusters);
    glutSetWindowTitle(title.c_str());
    
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    
    success = mesh.read(path);
    if (success) {
        generateClusterColors();
        mesh.segment(clusters);
    }
    
    printInstructions();
    glutInitWindowSize(gridX, gridY);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInit(&argc, argv);
    std::string title = "Shape Segmentation, clusters: " + std::to_string(clusters);
    glutCreateWindow(title.c_str());
    init();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMainLoop();
    
    return 0;
}
