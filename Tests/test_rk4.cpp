// test_rungekutta.cpp
#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "Solvers/SolverLibrary.hpp"
#include <cmath>
#include <vector>
#include "catch2/catch_all.hpp"
#include <iostream>

using Catch::Approx;

TEST_CASE("Runge-Kutta solver produces expected output for simple electromechanical system", "[solver][rungekutta]") {
    SolverLibrary solver;

    double t0 = 0.0;
    double tf = 1.0;
    double dt = 0.5;

    std::vector<double> y0 = {0.0, 1.0, 0.0, 0.0};  // initial: y, q, dy, dq

    std::vector<std::vector<double>> trajectory;
    int status = solver.solve_with_runge_kutta(y0, t0, tf, dt, trajectory);
    REQUIRE(status == 0);

    // Check size of trajectory
    REQUIRE(trajectory.size() == static_cast<size_t>((tf - t0) / dt) + 1);

    REQUIRE(Approx(trajectory[0][1]).margin(1e-2) == 0.0);         // y(0)
    REQUIRE(Approx(trajectory[1][1]).margin(1e-2) == 0.0552649927833715);    // y(0.5)
    REQUIRE(Approx(trajectory[2][1]).margin(1e-2) == 0.209457937508594);    // y(1.0)

    REQUIRE(Approx(trajectory[0][2]).margin(1e-2) == 1.0);         // q(0)
    REQUIRE(Approx(trajectory[1][2]).margin(1e-2) == 0.913098715067955);    // q(0.5)
    REQUIRE(Approx(trajectory[2][2]).margin(1e-2) == 0.660739021986070);    // q(1.0)

    // Optional: Check dy and dq as well
    REQUIRE(Approx(trajectory[1][3]).margin(1e-2) == 0.217138136315970);    // dy(0.5)
    REQUIRE(Approx(trajectory[2][3]).margin(1e-2) == 0.388411951826135);    // dy(1.0)

    REQUIRE(Approx(trajectory[1][4]).margin(1e-2) == -0.345108682051350);   // dq(0.5)
    REQUIRE(Approx(trajectory[2][4]).margin(1e-2) == -0.654109205664673);   // dq(1.0)
    
}  
