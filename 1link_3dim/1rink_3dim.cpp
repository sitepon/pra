#include <iostream>
#include <fstream>
#include <cmath>

#include <GL/glut.h>

struct Parameters
{
    // 重力加速度
    double g;

    // リンクの質量
    double m1;

    // 質量中心基準の慣性モーメント
    double I11;
    double I22;
    double I33;

};

// ===================
// グローバル変数
// ===================

Parameters p;

double y[12];

double dt = 0.0001;
double t_now = 0.0;
double t_end = 5.0;

// 剛体サイズ
double body_x = 1.0;
double body_y = 0.5;
double body_z = 5.0;

// ===================
// 運動方程式
// ===================
void func(double y[], double dydt[], const Parameters& p)
{
    // 一般化力
    double F1 = 0.0;
    double F2 = 0.0;
    double F3 = 0.0;

    double tau1 = 0.0; 
    double tau2 = 0.0; 
    double tau3 = 0.0; 

    // F_
    double F1_ = F1/p.m1;
    double F2_ = F2/p.m1;
    double F3_ = F3/p.m1 - p.g;

    double tau1_ = (tau1 - (p.I33 - p.I22)*y[10]*y[11])/p.I11; 
    double tau2_ = (tau2 - (p.I11 - p.I33)*y[9]*y[11])/p.I22; 
    double tau3_ = (tau3 - (p.I22 - p.I11)*y[9]*y[10])/p.I33;

    dydt[0] = y[6];
    dydt[1] = y[7];
    dydt[2] = y[8];
    dydt[3] = y[9] + sin(y[3])*tan(y[4])*y[10] + cos(y[3])*tan(y[4])*y[11];
    dydt[4] = cos(y[3])*y[10] - sin(y[3])*y[11];
    dydt[5] = sin(y[3])*y[10]/cos(y[4]) + cos(y[3])*y[11]/cos(y[4]);
    dydt[6] = F1_;
    dydt[7] = F2_;
    dydt[8] = F3_;
    dydt[9] = tau1_;
    dydt[10] = tau2_;
    dydt[11] = tau3_;
}

// ===================
// RK4
// ===================
void simulate()
{
    double k1[12];
    double k2[12];
    double k3[12];
    double k4[12];

    double y_tmp[12];

      // ルンゲクッタ法
        // k1
        func(y, k1, p);

        // k2
        y_tmp[0] = y[0] + 0.5 * dt * k1[0];
        y_tmp[1] = y[1] + 0.5 * dt * k1[1];
        y_tmp[2] = y[2] + 0.5 * dt * k1[2];
        y_tmp[3] = y[3] + 0.5 * dt * k1[3];
        y_tmp[4] = y[4] + 0.5 * dt * k1[4];
        y_tmp[5] = y[5] + 0.5 * dt * k1[5];
        y_tmp[6] = y[6] + 0.5 * dt * k1[6];
        y_tmp[7] = y[7] + 0.5 * dt * k1[7];
        y_tmp[8] = y[8] + 0.5 * dt * k1[8];
        y_tmp[9] = y[9] + 0.5 * dt * k1[9];
        y_tmp[10] = y[10] + 0.5 * dt * k1[10];
        y_tmp[11] = y[11] + 0.5 * dt * k1[11];

        func(y_tmp, k2, p);

        // k3
        y_tmp[0] = y[0] + 0.5 * dt * k2[0];
        y_tmp[1] = y[1] + 0.5 * dt * k2[1];
        y_tmp[2] = y[2] + 0.5 * dt * k2[2];
        y_tmp[3] = y[3] + 0.5 * dt * k2[3];
        y_tmp[4] = y[4] + 0.5 * dt * k2[4];
        y_tmp[5] = y[5] + 0.5 * dt * k2[5];
        y_tmp[6] = y[6] + 0.5 * dt * k2[6];
        y_tmp[7] = y[7] + 0.5 * dt * k2[7];
        y_tmp[8] = y[8] + 0.5 * dt * k2[8];
        y_tmp[9] = y[9] + 0.5 * dt * k2[9];
        y_tmp[10] = y[10] + 0.5 * dt * k2[10];
        y_tmp[11] = y[11] + 0.5 * dt * k2[11];

        func(y_tmp, k3, p);

        // k4
        y_tmp[0] = y[0] + dt * k3[0];
        y_tmp[1] = y[1] + dt * k3[1];
        y_tmp[2] = y[2] + dt * k3[2];
        y_tmp[3] = y[3] + dt * k3[3];
        y_tmp[4] = y[4] + dt * k3[4];
        y_tmp[5] = y[5] + dt * k3[5];
        y_tmp[6] = y[6] + dt * k3[6];
        y_tmp[7] = y[7] + dt * k3[7];
        y_tmp[8] = y[8] + dt * k3[8];
        y_tmp[9] = y[9] + dt * k3[9];
        y_tmp[10] = y[10] + dt * k3[10];
        y_tmp[11] = y[11] + dt * k3[11];

        func(y_tmp, k4, p);

        // 更新
        y[0] = y[0] + (dt/6) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);       // x
        y[1] = y[1] + (dt/6) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);       // y
        y[2] = y[2] + (dt/6) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);       // z
        y[3] = y[3] + (dt/6) * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]);       // φ
        y[4] = y[4] + (dt/6) * (k1[4] + 2.0 * k2[4] + 2.0 * k3[4] + k4[4]);       // θ
        y[5] = y[5] + (dt/6) * (k1[5] + 2.0 * k2[5] + 2.0 * k3[5] + k4[5]);       // Ψ
        y[6] = y[6] + (dt/6) * (k1[6] + 2.0 * k2[6] + 2.0 * k3[6] + k4[6]);       // x'
        y[7] = y[7] + (dt/6) * (k1[7] + 2.0 * k2[7] + 2.0 * k3[7] + k4[7]);       // y'
        y[8] = y[8] + (dt/6) * (k1[8] + 2.0 * k2[8] + 2.0 * k3[8] + k4[8]);       // z'
        y[9] = y[9] + (dt/6) * (k1[9] + 2.0 * k2[9] + 2.0 * k3[9] + k4[9]);       // w101
        y[10] = y[10] + (dt/6) * (k1[10] + 2.0 * k2[10] + 2.0 * k3[10] + k4[10]); // w102
        y[11] = y[11] + (dt/6) * (k1[11] + 2.0 * k2[11] + 2.0 * k3[11] + k4[11]); // w103
    
    t_now += dt;
}

// ===================
// 回転行列A10
// ===================
void applyRotation()
{
    double phi   = y[3];
    double theta = y[4];
    double psi   = y[5];

    double c1 = cos(phi);
    double s1 = sin(phi);

    double c2 = cos(theta);
    double s2 = sin(theta);

    double c3 = cos(psi);
    double s3 = sin(psi);


    /*
       Z-Y-X系

       R = Rz(psi) Ry(theta) Rx(phi)

    */

    double R[16] =
    {
        c3*c2,
        s3*c2,
        -s2,
        0.0,

        c3*s2*s1 - s3*c1,
        s3*s2*s1 + c3*c1,
        c2*s1,
        0.0,

        c3*s2*c1 + s3*s1,
        s3*s2*c1 - c3*s1,
        c2*c1,
        0.0,

        0.0,
        0.0,
        0.0,
        1.0
    };

    glMultMatrixd(R);
}

// ===================
// 座標軸
// ===================
void drawAxis()
{
    glLineWidth(2.0);

    glBegin(GL_LINES);

    // X軸
    glColor3f(1.0, 0.0, 0.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(2.0, 0.0, 0.0);

    // Y軸
    glColor3f(0.0, 1.0, 0.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 2.0, 0.0);

    // Z軸
    glColor3f(0.0, 0.0, 1.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 2.0);

    glEnd();
}

// ===================
// 剛体描画
// ===================
void drawBody()
{
    glPushMatrix();

    // ---------------------------
    // 重心位置へ並進
    // ---------------------------
    glTranslated(y[0], y[1], y[2]);

    // ---------------------------
    // 姿勢回転
    // ---------------------------
    applyRotation();

    // ---------------------------
    // 剛体座標系
    // ---------------------------
    glBegin(GL_LINES);
    // body x
    glColor3f(1.0, 0.0, 0.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);
    // body y
    glColor3f(0.0, 1.0, 0.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);

    // body z
    glColor3f(0.0, 0.0, 1.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);

    glEnd();

    // ---------------------------
    // 直方体
    // ---------------------------

    glColor3f(0.3, 0.6, 0.9);

    glScaled(body_x, body_y, body_z);

    glutSolidCube(1.0);

    glPopMatrix();
}

// ===================
// 地面
// ===================
void drawGround()
{
    glColor3f(0.7, 0.7, 0.7);

    glBegin(GL_LINES);

    for (double x = -5.0; x <= 5.0; x += 0.5)
    {
        glVertex3d(x, -5.0, 0.0);
        glVertex3d(x,  5.0, 0.0);
    }

    for (double y = -5.0; y <= 5.0; y += 0.5)
    {
        glVertex3d(-5.0, y, 0.0);
        glVertex3d( 5.0, y, 0.0);
    }

    glEnd();
}

// ===================
// 描画
// ===================
void display()
{
    // 前の画像と奥行き情報を全部消して、新しい1フレームを描く準備
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();

    // カメラ
    gluLookAt(
        5.0, -7.0, 4.0,   // カメラ位置

        0.0,  0.0, 0.0,   // 注視点

        0.0,  0.0, 1.0    // 上方向
    );

    drawGround();

    drawAxis();

    drawBody();

    glutSwapBuffers();
}

// ===================
// ウィンドウサイズ変更
// ===================
void reshape(int w, int h)
{
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();

    gluPerspective(
        45.0,                       // 視野角
        (double)w / (double)h,      // アスペクト比
        0.1,                        // 手前
        100.0                       // 奥
    );

    glMatrixMode(GL_MODELVIEW);
}

// ===================
// タイマー
// ===================
void timer(int value)
{
    int steps = 160;

    for (int i = 0; i < steps; i++)
    {
        simulate();
    }

    // 画面を更新してくださいと依頼
    glutPostRedisplay();

    glutTimerFunc(16, timer, 0);
}

int main(int argc, char** argv)
{
    // パラメータ
    p.g = 9.81;
    p.m1 = 1.0;
    p.I11 = 1.0;
    p.I22 = 1.0;
    p.I33 = 1.0;

    // 初期条件
    y[0] = 0.0;
    y[1] = 0.0;
    y[2] = 2.0;
    y[3] = 0.0;
    y[4] = 0.0;
    y[5] = 0.0;
    y[6] = 0.0;
    y[7] = 0.0;
    y[8] = 0.0;
    y[9] = 1.0;
    y[10] = 0.5;
    y[11] = 0.2;

    glutInit(&argc, argv);

    glutInitDisplayMode(
        GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH
    );

    glutInitWindowSize(
        800, 800
    );

    glutCreateWindow(
        "3D Rigid Body"
    );

    // ---------------------------
    // OpenGL設定
    // ---------------------------

    glClearColor(
        1.0, 1.0, 1.0, 1.0
    );

    glEnable(
        GL_DEPTH_TEST
    );

    // ---------------------------
    // コールバック
    // ---------------------------
    glutDisplayFunc(display);

    glutReshapeFunc(reshape);

    // 描画
    glutTimerFunc(16, timer, 0);

    glutMainLoop();

    return 0;

}