#include <iostream>
#include <vector>
#include "ODEs.cpp"
#include "SolverLibrary.hpp"


int SolverLibrary::solve_with_runge_kutta(std::vector<double>& y_initial, double t_start, double t_end, double t_step,
    std::vector<std::vector<double>>& trajectory) {
    const int N = 4;
    const double C1 = 0.2;
    const double C2 = 0.7;
    const double C3 = 1.0 - C2;

    std::vector<double> state = y_initial;
    std::vector<double> dydt(N);
    std::vector<double> k1(N), k2(N), k3(N), k4(N);

    double t = t_start;

    auto compute_derivatives = [&](const std::vector<double>& s, std::vector<double>& dsdt, double /*t*/) {
        dsdt[0] = s[2];                                          // dy/dt
        dsdt[1] = s[3];                                          // dq/dt
        dsdt[2] = 0.5 - s[0] - (C1 / 4.0) * std::pow(s[1], 2);   // d²y/dt²
        dsdt[3] = -C2 * s[1] - 2.0 * C3 * s[1] * s[0];           // d²q/dt²
    };

    trajectory.clear();
    trajectory.push_back({t, state[0], state[1], state[2], state[3]});  // Initial state

    while (t < t_end) {
        compute_derivatives(state, dydt, t);
        for (int i = 0; i < N; ++i) k1[i] = dydt[i] * t_step;

        std::vector<double> temp(N);
        for (int i = 0; i < N; ++i) temp[i] = state[i] + 0.5 * k1[i];
             compute_derivatives(temp, dydt, t + 0.5 * t_step);
        for (int i = 0; i < N; ++i) k2[i] = dydt[i] * t_step;

        for (int i = 0; i < N; ++i) temp[i] = state[i] + 0.5 * k2[i];
             compute_derivatives(temp, dydt, t + 0.5 * t_step);
        for (int i = 0; i < N; ++i) k3[i] = dydt[i] * t_step;

        for (int i = 0; i < N; ++i) temp[i] = state[i] + k3[i];
            compute_derivatives(temp, dydt, t + t_step);
        for (int i = 0; i < N; ++i) k4[i] = dydt[i] * t_step;

        for (int i = 0; i < N; ++i)
            state[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;

        t += t_step;
        trajectory.push_back({t, state[0], state[1], state[2], state[3]});
    }

    y_initial = state;
    return 0;
}
