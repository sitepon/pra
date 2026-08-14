#include <GL/freeglut.h>

void display()
{
    // 背景色
    glClearColor(1.0, 1.0, 1.0, 1.0);

    // 画面を消去
    glClear(GL_COLOR_BUFFER_BIT);

    // 線の太さ
    glLineWidth(5.0);

    // 黒
    glColor3f(0.0, 0.0, 0.0);

    // 線を描く
    glBegin(GL_LINES);

    glVertex2f(-0.5, 0.0);
    glVertex2f( 0.5, 0.0);

    glEnd();

    glutSwapBuffers();
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

    glutInitWindowSize(800, 600);

    glutCreateWindow("OpenGL Test");

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}