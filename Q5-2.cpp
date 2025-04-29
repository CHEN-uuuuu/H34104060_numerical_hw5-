#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// 微分方程組的右邊函數
vector<double> f(double t, const vector<double>& u) {
    double u1 = u[0];
    double u2 = u[1];
    double du1 = 9*u1 + 24*u2 + 5*cos(t) - (1.0/3)*sin(t);
    double du2 = -24*u1 - 52*u2 - 9*cos(t) + (1.0/3)*sin(t);
    return {du1, du2};
}

// 解析解
vector<double> exact_solution(double t) {
    double u1 = 2*exp(-3*t) - exp(-39*t) + (1.0/3)*cos(t);
    double u2 = -exp(-3*t) + 2*exp(-39*t) - (1.0/3)*cos(t);
    return {u1, u2};
}

// Runge-Kutta 4方法
void runge_kutta_4(double t0, double t_end, double h_init, const vector<double>& u0) {
    double h = h_init;  // 實際步長可能會調整
    vector<double> t_values;
    vector<vector<double>> u_values;

    double t = t0;
    vector<double> u = u0;

    t_values.push_back(t);
    u_values.push_back(u);

    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;  // 避免超過範圍
        }

        vector<double> k1 = f(t, u);
        vector<double> k2 = f(t + h/2, {u[0] + h/2 * k1[0], u[1] + h/2 * k1[1]});
        vector<double> k3 = f(t + h/2, {u[0] + h/2 * k2[0], u[1] + h/2 * k2[1]});
        vector<double> k4 = f(t + h, {u[0] + h * k3[0], u[1] + h * k3[1]});

        u[0] += (h/6)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        u[1] += (h/6)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        t += h;

        t_values.push_back(t);
        u_values.push_back(u);
    }

    // 顯示結果
    cout << fixed << setprecision(8);
    cout << "\n===== RK4 結果 (h = " << h_init << ") =====" << endl;
    cout << setw(5) << "t" << " | "
         << setw(12) << "u1_RK4" << " | "
         << setw(12) << "u1_exact" << " | "
         << setw(10) << "err_u1" << " || "
         << setw(12) << "u2_RK4" << " | "
         << setw(12) << "u2_exact" << " | "
         << setw(10) << "err_u2" << endl;
    cout << string(90, '-') << endl;

    for (size_t i = 0; i < t_values.size(); ++i) {
        vector<double> exact = exact_solution(t_values[i]);
        double err_u1 = fabs(u_values[i][0] - exact[0]);
        double err_u2 = fabs(u_values[i][1] - exact[1]);

        cout << setw(5) << setprecision(2) << fixed << t_values[i] << " | "
             << setw(12) << setprecision(8) << u_values[i][0] << " | "
             << setw(12) << exact[0] << " | "
             << setw(10) << scientific << err_u1 << " || "
             << setw(12) << fixed << u_values[i][1] << " | "
             << setw(12) << exact[1] << " | "
             << setw(10) << scientific << err_u2 << endl;
    }
}

int main() {
    vector<double> u0 = {4.0/3.0, 2.0/3.0};
    double t0 = 0.0;
    double t_end = 1.0;

    // 執行不同步長
    runge_kutta_4(t0, t_end, 0.1, u0);
    runge_kutta_4(t0, t_end, 0.05, u0);

    return 0;
}
