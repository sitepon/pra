#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <GL/freeglut.h>

struct Parameters
{
    // 重力加速度
    double g;

    // 各リンクの質量
    double m1;
    double m2;

    // 各リンクの質量中心基準の慣性モーメント
    double I1;
    double I2;

    double l1;
    double l2;

    double r1;
    double r2;
};

struct State
{
    double t;
    double th1;
    double th2;
    double dth1;
    double dth2;
};

// グローバル変数

Parameters p;

std::vector<State> result;

// 現在表示しているフレーム
int frame = 0;

// 描画用の時間間隔
const double animation_dt = 1.0/60.0;

void func(double y[], double dydt[], const Parameters& p)
{
    // 質量行列 M
    double M11 = p.I1 + p.I2 + p.m1*p.r1*p.r1 + p.m2*(p.l1*p.l1 + p.r2*p.r2 + 2*p.l1*p.r2*cos(y[1]));
    double M12 = p.I2 + p.m2*(p.r2*p.r2 + p.l1*p.r2*cos(y[1]));
    double M21 = p.I2 + p.m2*(p.r2*p.r2 + p.l1*p.r2*cos(y[1]));
    double M22 = p.I2 + p.m2*p.r2*p.r2;

    double detM = M11*M22 - M12*M21;

    // 慣性力項
    double C1 = -p.m2*p.l1*p.r2*sin(y[1])*(2*y[2]*y[3] + y[3]*y[3]);
    double C2 = p.m2*p.l1*p.r2*sin(y[1])*y[2]*y[2];

    // 重力項
    double g1 = (p.m1*p.r1 + p.m2*p.l1)*p.g*cos(y[0]) + p.m2*p.g*p.r2*cos(y[0] + y[1]);
    double g2 = p.m2*p.g*p.r2*cos(y[0] + y[1]);

    // 一般化力
    double F1 = 0;
    double F2 = 0.0;

    // F_
    double F1_ = F1 - C1 - g1;
    double F2_ = F2 - C2 - g2;

    dydt[0] = y[2];
    dydt[1] = y[3];
    dydt[2] = (M22*F1_ - M12*F2_)/detM;
    dydt[3] = (-M21*F1_ + M11*F2_)/detM;
}

void simulate()
{
    // 初期条件
    double y[4];
    
    y[0] = 0.0;
    y[1] = 0.0;
    y[2] = 0.0;
    y[3] = 0.0;

    // 時間刻み
    double dt = 0.0001;
    double t_end = 5.0;

    // RK4
    double k1[4];
    double k2[4];
    double k3[4];
    double k4[4];

    double y_tmp[4];

    // 結果をcsvに保存
    std::ofstream file("result_2rink.csv");
    file << "t,x,z,th1,th2,dth1,dth2\n";

    double next_save_time = 0.0;

    for (double t = 0.0; t <= t_end; t += dt)
    {

        // --------------------------------
        // 60fps相当でアニメーション用データ保存
        // --------------------------------

        if (t >= next_save_time)
        {
            result.push_back(
            {
                t,
                y[0],
                y[1],
                y[2],
                y[3]
            });

            next_save_time += animation_dt;
        }

        double x = p.l1*cos(y[0]) + p.l2*cos(y[0] + y[1]);
        double z = p.l1*sin(y[0]) + p.l2*sin(y[0] + y[1]);

        // 結果保存
        file << t << ","
             << x << "," 
             << z << "," 
             << y[0] << "," 
             << y[1] << "," 
             << y[2] << "," 
             << y[3] << "\n";
        
        // ルンゲクッタ法

        // k1
        func(y, k1, p);

        // k2
        y_tmp[0] = y[0] + 0.5 * dt * k1[0];
        y_tmp[1] = y[1] + 0.5 * dt * k1[1];
        y_tmp[2] = y[2] + 0.5 * dt * k1[2];
        y_tmp[3] = y[3] + 0.5 * dt * k1[3];

        func(y_tmp, k2, p);

        // k3
        y_tmp[0] = y[0] + 0.5 * dt * k2[0];
        y_tmp[1] = y[1] + 0.5 * dt * k2[1];
        y_tmp[2] = y[2] + 0.5 * dt * k2[2];
        y_tmp[3] = y[3] + 0.5 * dt * k2[3];

        func(y_tmp, k3, p);

        // k4
        y_tmp[0] = y[0] + dt * k3[0];
        y_tmp[1] = y[1] + dt * k3[1];
        y_tmp[2] = y[2] + dt * k3[2];
        y_tmp[3] = y[3] + dt * k3[3];

        func(y_tmp, k4, p);

        // 更新
        y[0] = y[0] + (dt/6) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
        y[1] = y[1] + (dt/6) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
        y[2] = y[2] + (dt/6) * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);
        y[3] = y[3] + (dt/6) * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]);
    }

    file.close();

    std::cout << "計算終了" << std::endl;
    std::cout << "アニメーションフレーム数 = "
              << result.size()
              << std::endl;
}

void drawCircle(double x, double z, double radius)
{
    const int N = 30;

    glBegin(GL_POLYGON);

    for (int i = 0; i < N; i++)
    {
        double theta =
            2.0 * 3.1415 * i / N;

        glVertex2d(
            x + radius * cos(theta),
            z + radius * sin(theta)
        );
    }

    glEnd();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);

    if (result.empty())
    {
        glutSwapBuffers();
        return;
    }

    State s = result[frame];

    double th1 = s.th1;
    double th2 = s.th2;

    // ------------------------------
    // 順運動学
    // ------------------------------

    // 第1関節
    double x0 = 0.0;
    double z0 = 0.0;

    // 第2関節
    double x1 = p.l1*cos(th1);
    double z1 = p.l1*sin(th1);
    
    // 先端
    double x2 = x1 + p.l2 * cos(th1 + th2);
    double z2 = z1 + p.l2 * sin(th1 + th2);

    // ==============================
    // リンク描画
    // ==============================

    glLineWidth(8.0);

    glColor3f(0.1, 0.1, 0.1);

    glBegin(GL_LINES);

    // Link 1
    glVertex2d(x0, z0);
    glVertex2d(x1, z1);

    // Link 2
    glVertex2d(x1, z1);
    glVertex2d(x2, z2);

    glEnd();

    // ==============================
    // 関節描画
    // ==============================

    // 第1関節
    glColor3f(0.2, 0.2, 0.8);

    drawCircle(x0, z0, 0.07);

    // 第2関節
    glColor3f(0.8, 0.2, 0.2);

    drawCircle(x1, z1, 0.07);

    // 先端
    glColor3f(0.2, 0.7, 0.2);

    drawCircle(x2, z2, 0.07);

    glutSwapBuffers();
}

void timer(int value)
{
    frame++;


    // 最後まで行ったら最初に戻る
    if (frame >= static_cast<int>(result.size()))
    {
        frame = 0;
    }


    // 再描画を要求
    glutPostRedisplay();


    // 約60fps
    glutTimerFunc(
        16,
        timer,
        0
    );
}

void reshape(int width, int height)
{
    glViewport(
        0,
        0,
        width,
        height
    );


    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();


    // 描画範囲
    double range = 2.5;


    if (width >= height)
    {
        double aspect =
            static_cast<double>(width)
            / static_cast<double>(height);

        glOrtho(
            -range * aspect,
             range * aspect,
            -range,
             range,
            -1.0,
             1.0
        );
    }
    else
    {
        double aspect =
            static_cast<double>(height)
            / static_cast<double>(width);

        glOrtho(
            -range,
             range,
            -range * aspect,
             range * aspect,
            -1.0,
             1.0
        );
    }

    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
}

int main(int argc, char** argv)
{
    // =======================
    // パラメータ
    // =======================
    p.g = 9.81;

    p.m1 = 1.0;
    p.m2 = 1.0;

    p.I1 = 1.0;
    p.I2 = 1.0;

    p.l1 = 1.0;
    p.l2 = 1.0;

    p.r1 = 0.5;
    p.r2 = 0.5;

    // =======================
    // 数値計算
    // =======================
    simulate();

    // =======================
    // OpenGL 初期化
    // =======================
    glutInit(
        &argc,argv
    );

    glutInitDisplayMode(
        GLUT_DOUBLE | GLUT_RGB
    );
    
    glutInitWindowSize(
        800,800
    );

    glutCreateWindow(
        "2 Link Robot Animation"
    );

    // 背景色
    glClearColor(
        1.0,
        1.0,
        1.0,
        1.0
    );

    // コールバック関数登録
    glutDisplayFunc(display);

    glutReshapeFunc(reshape);

    glutTimerFunc(
        16,timer,0
    );

    // OpenGLメインループ
    glutMainLoop();

    return 0;
}