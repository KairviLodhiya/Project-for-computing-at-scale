#include <iostream>
#include <boost/numeric/odeint.hpp>
#include "SolverLibrary.hpp"

int SolverLibrary::solve_with_boost_odeint(std::vector<double>& y_initial, double t_start, double t_end, double t_step) {
    using namespace boost::numeric::odeint;

    std::vector<double> state = y_initial;
    runge_kutta4<std::vector<double>> stepper;

    double t = t_start;
    while (t < t_end) {
        stepper.do_step(SolverLibrary::derivatives, state, t, t_step);
        t += t_step;
    }

    return 0;
}