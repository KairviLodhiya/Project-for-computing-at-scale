#include "SolverLibrary.hpp"
#include <iostream>
#include <vector>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode.h>


void SolverLibrary::derivatives(const std::vector<double>& state, std::vector<double>& dydt, double t) {
    const double C1 = 0.2;
    const double C2 = 0.7;
    const double C3 = 1 - C2;

    double y1 = state[0];
    double y2 = state[1];
    double q1 = state[2];
    double q2 = state[3];

    dydt[0] = y2;
    dydt[1] = (0.5 - y1 - (C1 / 4.0) * q1 * q1);
    dydt[2] = q2;
    dydt[3] = -C2 * q1 - 2 * C3 * q1 * y1;
}

// int cv_derivatives(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
//     SolverLibrary* solver = static_cast<SolverLibrary*>(user_data);  // Access the SolverLibrary instance
//     std::vector<double> state(NV_LENGTH_S(y));  // Copy N_Vector to std::vector
//     for (int i = 0; i < NV_LENGTH_S(y); ++i) {
//         state[i] = NV_Ith_S(y, i);
//     }

//     std::vector<double> dydt(state.size());

//     solver->SolverLibrary::derivatives(state, dydt, t);  // Call the actual member function

//     // Copy the results back to N_Vector
//     for (int i = 0; i < dydt.size(); ++i) {
//         NV_Ith_S(ydot, i) = dydt[i];
//     }

//     return 0;  // Return 0 to indicate successful computation
// }


