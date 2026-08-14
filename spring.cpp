#include <iostream>
#include <fstream>

void func(double y[], double dydt[], double m, double k)
{
    dydt[0] = y[1];
    dydt[1] = -(k / m) * y[0];
}

int main()
{
    // パラメータ
    double m = 1.0;
    double k = 10.0;

    // 初期条件
    double y[2];
    
    y[0] = 0.1;
    y[1] = 0.0;

    // 時間刻み
    double dt = 0.001;
    double t_end = 10.0;

    // RK4
    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];

    double y_tmp[2];

    // 結果をcsvに保存
    std::ofstream file("result.csv");
    file << "t,x,v\n";

    for (double t = 0.0; t <= t_end; t += dt)
    {
        // 結果保存
        file << t << "," << y[0] << "," << y[1] << "\n";
        
        // ルンゲクッタ法

        // k1
        func(y, k1, m, k);

        // k2
        y_tmp[0] = y[0] + 0.5 * dt * k1[0];
        y_tmp[1] = y[1] + 0.5 * dt * k1[1];

        func(y_tmp, k2, m, k);

        // k3
        y_tmp[0] = y[0] + 0.5 * dt * k2[0];
        y_tmp[1] = y[1] + 0.5 * dt * k2[1];

        func(y_tmp, k3, m, k);

        // k4
        y_tmp[0] = y[0] + dt * k3[0];
        y_tmp[1] = y[1] + dt * k3[1];

        func(y_tmp, k4, m, k);

        // 更新
        y[0] = y[0] + (dt/6) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
        y[1] = y[1] + (dt/6) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
    }

    file.close();

    std::cout << "計算終了" << std::endl;

    return 0;

}