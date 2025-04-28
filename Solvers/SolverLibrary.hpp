// SolverLibrary.hpp
#ifndef SOLVERLIBRARY_HPP
#define SOLVERLIBRARY_HPP

#include <vector>
#include <sundials/include/cvode/cvode.h>                 // for CVode, N_Vector
#include <sundials/include/nvector/nvector_serial.h>      // for N_VSerial types
#include <sundials/sundials_config.h>

class SolverLibrary {
public:
    // Constructor
    SolverLibrary();

    // Solvers
    int solve_with_sundials(std::vector<double>& y_initial, double t_start, double t_end, double t_step, std::vector<std::vector<double>>& trajectory);
    int solve_with_boost_odeint(std::vector<double>& y_initial, double t_start, double t_end, double t_step);
    int solve_with_runge_kutta(std::vector<double>& y_initial, double t_start, double t_end, double t_step, std::vector<std::vector<double>>& trajectory);
    // Derivatives function
    static void derivatives(const std::vector<double>& state, std::vector<double>& dydt, double t);
    // private:
    // Static wrapper for Sundials function signature
     static int cv_derivatives(double t, N_Vector y, N_Vector ydot, void *user_data);
};

#endif // SOLVERLIBRARY_HPP
