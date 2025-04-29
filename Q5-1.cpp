#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// 定義微分方程
double f(double t, double y) {
    return 1 + (y / t) + pow(y / t, 2);
}

// 真實解
double exact_solution(double t) {
    return t * tan(log(t));
}

// 計算 f_t
double f_t(double t, double y) {
    double fy = (1.0 / t) + (2.0 * y / (t * t));
    double ft = -(y / (t * t)) - (2.0 * y * y / pow(t, 3));
    return ft + fy * f(t, y);
}

int main() {
    // 初始條件
    double t0 = 1.0;
    double y0 = 0.0;
    double h = 0.1;
    double t_end = 2.0;

    // 建立時間點
    int n = static_cast<int>((t_end - t0) / h) + 1;
    vector<double> t_values(n);
    for (int i = 0; i < n; ++i) {
        t_values[i] = t0 + i * h;
    }

    // Euler方法近似解
    vector<double> y_euler(n, 0.0);
    y_euler[0] = y0;

    for (int i = 1; i < n; ++i) {
        y_euler[i] = y_euler[i-1] + h * f(t_values[i-1], y_euler[i-1]);
    }

    // Taylor 2nd order 方法近似解
    vector<double> y_taylor2(n, 0.0);
    y_taylor2[0] = y0;

    for (int i = 1; i < n; ++i) {
        y_taylor2[i] = y_taylor2[i-1] 
                     + h * f(t_values[i-1], y_taylor2[i-1]) 
                     + (h * h / 2.0) * f_t(t_values[i-1], y_taylor2[i-1]);
    }

    // 真實解
    vector<double> y_exact(n);
    for (int i = 0; i < n; ++i) {
        y_exact[i] = exact_solution(t_values[i]);
    }

    // 顯示結果
    cout << fixed << setprecision(6);

    cout << "t\t\tEuler\t\tExact\t\tError (Euler)\n";
    for (int i = 0; i < n; ++i) {
        cout << t_values[i] << "\t" << y_euler[i] << "\t" << y_exact[i] << "\t" << fabs(y_exact[i] - y_euler[i]) << "\n";
    }

    cout << "\n------------------------------------------------------------\n" << endl;

    cout << "t\t\tTaylor 2nd\tExact\t\tError (Taylor 2nd)\n";
    for (int i = 0; i < n; ++i) {
        cout << t_values[i] << "\t" << y_taylor2[i] << "\t" << y_exact[i] << "\t" << fabs(y_exact[i] - y_taylor2[i]) << "\n";
    }

    return 0;
}
