#include "SolverLibrary.hpp"
#include <iostream>

#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_context.h>

int SolverLibrary::cv_derivatives(double t, N_Vector y, N_Vector ydot, void* user_data) {
    double* Y = NV_DATA_S(y);
    double* dY = NV_DATA_S(ydot);

    const double C1 = 0.2;
    const double C2 = 0.7;
    const double C3 = 1.0 - C2;

    dY[0] = Y[2];
    dY[1] = Y[3];
    dY[2] = 0.5 - Y[0] -( (C1 / 4.0) * pow(Y[1],2));
    dY[3] = -(C2 * Y[1] )- (2.0 * C3 * Y[1] * Y[0]);

    return 0;
}

SolverLibrary::SolverLibrary() {}

int SolverLibrary::solve_with_sundials(std::vector<double>& y_initial, double t0, double tf, double dt, std::vector<std::vector<double>>& trajectory)
{
    // Step 1: Create context
    SUNContext sunctx;
    if (SUNContext_Create(NULL, &sunctx)) {
        std::cerr << "[ERROR] SUNContext creation failed.\n";
        return -1;
    }

    // Step 2: Create vector
    N_Vector y = N_VNew_Serial(4, sunctx);
    if (!y) {
        std::cerr << "[ERROR] N_Vector allocation failed.\n";
        return -1;
    }
    for (int i = 0; i < 4; ++i){
        NV_Ith_S(y, i) = y_initial[i];
        
    }

    // Step 3: Create CVODE memory
    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (!cvode_mem) {
        std::cerr << "[ERROR] CVodeCreate failed.\n";
        return -1;
    }

    // Step 4: Initialize CVODE and set tolerances
    if (CVodeInit(cvode_mem, cv_derivatives, t0, y) != CV_SUCCESS) {
        std::cerr << "[ERROR] CVodeInit failed.\n";
        return -1;
    }

    if (CVodeSStolerances(cvode_mem, 1e-6, 1e-7) != CV_SUCCESS) {
        std::cerr << "[ERROR] Tolerance setup failed.\n";
        return -1;
    }

    // Step 5: Optional integration controls
    CVodeSetInitStep(cvode_mem, 1e-10);    // Start small
    CVodeSetMaxStep(cvode_mem, 1);      // Adaptive range
    CVodeSetMaxNumSteps(cvode_mem, 500000000);

    // Step 6: Create SUNMatrix and Linear Solver
    SUNMatrix A = SUNDenseMatrix(4, 4, sunctx);
    if (!A) {
        std::cerr << "[ERROR] SUNDenseMatrix failed.\n";
        return -1;
    }

    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    if (!LS) {
        std::cerr << "[ERROR] SUNLinSol_Dense failed.\n";
        return -1;
    }

    if (CVodeSetLinearSolver(cvode_mem, LS, A) != CV_SUCCESS) {
        std::cerr << "[ERROR] CVodeSetLinearSolver failed.\n";
        return -1;
    }

    // Step 7: Nonlinear solver (Newton)
    SUNNonlinearSolver NLS = SUNNonlinSol_Newton(y, sunctx);
    if (!NLS) {
        std::cerr << "[ERROR] SUNNonlinSol_Newton failed.\n";
        return -1;
    }

    if (CVodeSetNonlinearSolver(cvode_mem, NLS) != CV_SUCCESS) {
        std::cerr << "[ERROR] CVodeSetNonlinearSolver failed.\n";
        return -1;
    }

    // double t = t0;
    // double tout = t0 + dt_unused;

    // while (t < tf) {
    //     if (tout > tf) tout = tf;
    
    //     int flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    //     if (flag < 0) {
    //         std::cerr << "CVode failed at t = " << t << std::endl;
    //         break;
    //     }
    
    //     tout += dt_unused;
    // }
    // for (int i = 0; i < 4; ++i) y_initial[i] = NV_Ith_S(y, i);

    // Step 8: Integration loop
    double t = t0;
    double t_out = t0 + dt;

    std::vector<double> state0(5);
    state0[0] = t0;
    state0[1] = NV_Ith_S(y, 0);
    state0[2] = NV_Ith_S(y, 1);
    state0[3] = NV_Ith_S(y, 2);
    state0[4] = NV_Ith_S(y, 3);
    trajectory.push_back(state0);

    while (t <= tf) {

        int flag = CVode(cvode_mem, t_out, y, &t, CV_NORMAL);
        if (flag < 0) {
            std::cerr << "[ERROR] CVode failed at t = " << t << " with flag = " << flag << "\n";
            break;
        }
        std::vector<double> state(5);
        state[0] = t;               // actual time
        state[1] = NV_Ith_S(y, 0);  // y
        state[2] = NV_Ith_S(y, 1);  // q
        state[3] = NV_Ith_S(y, 2);  // dy
        state[4] = NV_Ith_S(y, 3);  // dq
        // std::cout << "t = " << t << std::endl;
        //         << ", y = " << state[1]
        //         << ", q = " << state[2]
        //         << ", dy = " << state[3]
        //         << ", dq = " << state[4] << std::endl;
    
        trajectory.push_back(state);
        
        
        t_out += dt;
     }
    

    for (int i = 0; i < 4; ++i) y_initial[i] = NV_Ith_S(y, i);
    

    // Step 9: Clean-up
    CVodeFree(&cvode_mem);
    N_VDestroy(y);
    SUNMatDestroy(A);
    SUNLinSolFree(LS);
    SUNNonlinSolFree(NLS);
    SUNContext_Free(&sunctx);

    return 0;
}
