#define CATCH_CONFIG_MAIN
// #include <catch2/catch.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "Solvers/SolverLibrary.hpp"
#include "catch2/catch_all.hpp"
#include <iostream>

using Catch::Approx;

TEST_CASE("ODE solver returns expected shape", "[trajectory]") {
    SolverLibrary solver;
    std::vector<double> y0 = {0.0, 1.0, 0.0, 0.0}; // Initial conditions [y, q, dy, dq]
    double t0 = 0.0;
    double tf = 10.0;
    double dt = 0.001;

    std::vector<std::vector<double>> trajectory;
    int status = solver.solve_with_sundials(y0, t0, tf, dt, trajectory);

    REQUIRE(status == 0);
    REQUIRE(trajectory.size()-1 == static_cast<size_t>((tf - t0) / dt) + 1);
    REQUIRE(trajectory[0].size() == 5); // t, y, q, dy, dq
}

TEST_CASE("Check initial condition preservation", "[trajectory]") {
    SolverLibrary solver;
    std::vector<double> y0 = {0.5, 1.0, 0.0, 0.0};
    double t0 = 0.0;
    double tf = 1.0;
    double dt = 0.1;

    std::vector<std::vector<double>> trajectory;
    solver.solve_with_sundials(y0, t0, tf, dt, trajectory);

    REQUIRE(trajectory.front()[1] == Approx(0.5).margin(1e-10)); // y
    REQUIRE(trajectory.front()[2] == Approx(1.0).margin(1e-10)); // q
    REQUIRE(trajectory.front()[3] == Approx(0.0).margin(1e-10)); // dy
    REQUIRE(trajectory.front()[4] == Approx(0.0).margin(1e-10)); // dq
}

TEST_CASE("Check accuracy of results", "[trajectory]") {
    SolverLibrary solver;
    std::vector<double> y0 = {0, 1.0, 0.0, 0.0};

    double t0 = 0.0;
    double tf = 1.0;
    double dt = 0.5;

    std::vector<std::vector<double>> trajectory;
    solver.solve_with_sundials(y0, t0, tf, dt, trajectory);

    int status = solver.solve_with_sundials(y0, t0, tf, dt, trajectory);
    REQUIRE(status == 0);

    // Check size of trajectory
    //REQUIRE(trajectory.size()-1 == static_cast<size_t>((tf - t0) / dt) + 1);
    

    // Extract and test specific values at t = 0, 0.5, and 1.0
    REQUIRE(Approx(trajectory[0][1]).margin(1e-4) == 0.0);         // y(0)
    REQUIRE(Approx(trajectory[1][1]).margin(1e-4) == 0.0552649927833715);    // y(0.5)
    REQUIRE(Approx(trajectory[2][1]).margin(1e-4) == 0.209457937508594);    // y(1.0)

    REQUIRE(Approx(trajectory[0][2]).margin(1e-4) == 1.0);         // q(0)
    REQUIRE(Approx(trajectory[1][2]).margin(1e-4) == 0.913098715067955);    // q(0.5)
    REQUIRE(Approx(trajectory[2][2]).margin(1e-4) == 0.660739021986070);    // q(1.0)

    // Optional: Check dy and dq as well
    REQUIRE(Approx(trajectory[1][3]).margin(1e-4) == 0.217138136315970);    // dy(0.5)
    REQUIRE(Approx(trajectory[2][3]).margin(1e-4) == 0.388411951826135);    // dy(1.0)

    REQUIRE(Approx(trajectory[1][4]).margin(1e-4) == -0.345108682051350);   // dq(0.5)
    REQUIRE(Approx(trajectory[2][4]).margin(1e-4) == -0.654109205664673);   // dq(1.0)
}
