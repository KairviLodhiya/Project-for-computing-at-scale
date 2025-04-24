#include <iostream>
#include <vector>
#include "ODEs.cpp"
#include "SolverLibrary.hpp"


int SolverLibrary::solve_with_runge_kutta(std::vector<double>& y_initial, double t_start, double t_end, double t_step) {
    const int N = 4;
    std::vector<double> state = y_initial;
    std::vector<double> dydt(N);
    std::vector<double> k1(N), k2(N), k3(N), k4(N);

    double t = t_start;
    while (t < t_end) {
        derivatives(state, dydt, t);
        for (int i = 0; i < N; ++i) {
            k1[i] = dydt[i] * t_step;
        }

        std::vector<double> temp_state(N);
        for (int i = 0; i < N; ++i) {
            temp_state[i] = state[i] + 0.5 * k1[i];
        }
        derivatives(temp_state, dydt, t + 0.5 * t_step);
        for (int i = 0; i < N; ++i) {
            k2[i] = dydt[i] * t_step;
        }

        for (int i = 0; i < N; ++i) {
            temp_state[i] = state[i] + 0.5 * k2[i];
        }
        derivatives(temp_state, dydt, t + 0.5 * t_step);
        for (int i = 0; i < N; ++i) {
            k3[i] = dydt[i] * t_step;
        }

        for (int i = 0; i < N; ++i) {
            temp_state[i] = state[i] + k3[i];
        }
        derivatives(temp_state, dydt, t + t_step);
        for (int i = 0; i < N; ++i) {
            k4[i] = dydt[i] * t_step;
        }

        // Update the state using the Runge-Kutta method
        for (int i = 0; i < N; ++i) {
            state[i] += (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0;
        }

        t += t_step;
    }

    return 0;
}