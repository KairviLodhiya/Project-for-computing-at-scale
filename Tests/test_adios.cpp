#define CATCH_CONFIG_MAIN
// #include <catch2/catch.hpp>
#include <adios2.h>
#include <filesystem>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "catch2/catch_all.hpp"
#include "Solvers/SolverLibrary.hpp"

using Catch::Approx;

TEST_CASE("ADIOS2 Save and Read", "[adios2]") {
    std::string bp_file = "test_output.bp";

    // Step 1: Write some data
    {
        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("TestIO");
        adios2::Engine writer = io.Open(bp_file, adios2::Mode::Write);

        std::vector<double> time = {0.0, 0.1, 0.2};
        std::vector<double> y    = {0.0, 0.5, 1.0};
        std::vector<double> q    = {1.0, 0.8, 0.6};
        std::vector<double> dy   = {0.0, -0.2, -0.4};
        std::vector<double> dq   = {0.0,  0.1,  0.2};

        auto v_time = io.DefineVariable<double>("time", {}, {}, {time.size()});
        auto v_y    = io.DefineVariable<double>("y", {}, {}, {y.size()});
        auto v_q    = io.DefineVariable<double>("q", {}, {}, {q.size()});
        auto v_dy   = io.DefineVariable<double>("dy", {}, {}, {dy.size()});
        auto v_dq   = io.DefineVariable<double>("dq", {}, {}, {dq.size()});

        writer.Put(v_time, time.data());
        writer.Put(v_y,    y.data());
        writer.Put(v_q,    q.data());
        writer.Put(v_dy,   dy.data());
        writer.Put(v_dq,   dq.data());
        writer.Close();
    }

    // Step 2: Verify file exists
    REQUIRE(std::filesystem::exists(bp_file));

    // Step 3: Read back the data
    {
        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("ReaderIO");
        adios2::Engine reader = io.Open(bp_file, adios2::Mode::Read);

        REQUIRE(reader.BeginStep() == adios2::StepStatus::OK);

        auto vars = io.AvailableVariables();
        REQUIRE(vars.count("time") > 0);
        REQUIRE(vars.count("y") > 0);
        REQUIRE(vars.count("q") > 0);
        REQUIRE(vars.count("dy") > 0);
        REQUIRE(vars.count("dq") > 0);

        reader.EndStep();
        reader.Close();
    }

    // Cleanup
    if (std::filesystem::exists(bp_file)) {
        std::filesystem::remove_all(bp_file);
    }
}
