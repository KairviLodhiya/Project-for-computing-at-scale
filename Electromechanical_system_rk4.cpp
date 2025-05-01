
#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <adios2.h>
#include "Solvers/SolverLibrary.hpp"
#include <chrono>

// Struct that defines the ODE system
typedef double real;

struct StateSnapshot {
  real time;
  real y, dy, q, dq;
};

/*struct ElectromechanicalSystem {
  double C1, C2, C3;

  KOKKOS_INLINE_FUNCTION
  ElectromechanicalSystem(double C1_in, double C2_in)
    : C1(C1_in), C2(C2_in), C3(1.0 - C2_in) {}

  // y = [y, y_dot, q, q_dot]
  KOKKOS_INLINE_FUNCTION
  void operator()(const double* y, double* dydt) const {
    dydt[0] = y[1]; // dy/dt = y'
    dydt[1] = 0.5 - y[0] - (C1 / 4.0) * y[2] * y[2];
    dydt[2] = y[3]; // dq/dt = q'
    dydt[3] = -C2 * y[2] - 2.0 * C3 * y[2] * y[0];
  }
};*/

// RK4 time stepping
// KOKKOS_INLINE_FUNCTION
/* void rk4_step(const ElectromechanicalSystem& sys, double* y, double dt) {
  double dydt[4], k1[4], k2[4], k3[4], k4[4], y_temp[4];
  sys(y, k1);
  for (int i = 0; i < 4; ++i) y_temp[i] = y[i] + 0.5 * dt * k1[i];
  sys(y_temp, k2);
  for (int i = 0; i < 4; ++i) y_temp[i] = y[i] + 0.5 * dt * k2[i];
  sys(y_temp, k3);
  for (int i = 0; i < 4; ++i) y_temp[i] = y[i] + dt * k3[i];
  sys(y_temp, k4);
  for (int i = 0; i < 4; ++i)
    y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}*/

int main(int argc, char* argv[]) {

  auto start_time = std::chrono::high_resolution_clock::now();
    const int N = 10;
    const real t_start = 0.0;
    const real t_end = 5000;
    const real dt = 0.001;

    std::vector<std::vector<StateSnapshot>> all_simulations(N);

    SolverLibrary solver;

    adios2::ADIOS adios;
    std::string output_file = "all_groups_rk4.bp";
    adios2::IO io = adios.DeclareIO("simulation_io");
    adios2::Engine writer = io.Open(output_file, adios2::Mode::Write);
    

    for (int i = 0; i < N; ++i) {
      
      double y_init = i * 0.1;
      double q_init = 1.0;
      
      std::vector<double> y0 = {y_init, q_init, 0.0 , 0.0};
      
      std::vector<std::vector<double>> y_out;
      std::vector<std::vector<double>> trajectory;

      int status = solver.solve_with_runge_kutta(y0, t_start, t_end, dt, trajectory);
      
      std::vector<StateSnapshot> snapshots;
      std::cout << "[CHECK] trajectory size = " << trajectory.size() << ", expected = " << static_cast<int>((t_end - t_start) / dt) + 1 << std::endl;
      StateSnapshot snap;
      
      for (int t = 0;  t <= trajectory.size()-1; ++t) {
        
        
        snap.time = trajectory[t][0];
        snap.y  = trajectory[t][1];
        snap.q = trajectory[t][2];
        snap.dy  = trajectory[t][3];
        snap.dq = trajectory[t][4];
        snapshots.push_back(snap);

        // std::cout << "t = " << snap.time << std::endl;
      //           << ", y = " << snap.y
      //           << ", q = " << snap.q
      //           << ", dy = " << snap.dy
      //           << ", dq = " << snap.dq << std::endl;
       }
      
      int y_index = static_cast<int>(y_init * 10);
      int q_index = static_cast<int>(q_init * 100);
      
      std::string group_name = "y" + std::to_string(y_index) + "_q" + std::to_string(q_index);
      
      std::cout << "Group name = " << group_name << "\n";
      std::cout << "[checkpoint] Writing to ADIOS" << std::endl;
      

      std::vector<real> time_vec, y_vec, dy_vec, q_vec, dq_vec;
      for (const auto& row : trajectory) {
          time_vec.push_back(row[0]);  // t
          y_vec.push_back(row[1]);     // y
          q_vec.push_back(row[2]);     // q
          dy_vec.push_back(row[3]);    // dy
          dq_vec.push_back(row[4]);    // dq
    }
      
      auto v_time = io.DefineVariable<real>(group_name + "/time", {}, {}, {time_vec.size()});
      auto v_y    = io.DefineVariable<real>(group_name + "/y",    {}, {}, {y_vec.size()});
      auto v_q    = io.DefineVariable<real>(group_name + "/q",    {}, {}, {q_vec.size()});
      auto v_dy   = io.DefineVariable<real>(group_name + "/dy",   {}, {}, {dy_vec.size()});
      auto v_dq   = io.DefineVariable<real>(group_name + "/dq",   {}, {}, {dq_vec.size()});

      if (!v_time || !v_y || !v_q || !v_dy || !v_dq) {
        std::cerr << "[ERROR] Variable definition failed in group: " << group_name << std::endl;
        writer.Close();
        continue;
    }
      writer.BeginStep();
      writer.Put(v_time, time_vec.data());
      writer.Put(v_y,    y_vec.data());
      writer.Put(v_dy,   dy_vec.data());
      writer.Put(v_q,    q_vec.data());
      writer.Put(v_dq,   dq_vec.data());
      writer.EndStep();
      
      std::cout << "[checkpoint] Saved the file" << std::endl;
    }
    writer.Close();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << "[DONE] All simulations completed and saved." << std::endl;
    std::cout << "[TIME] Total execution time: " << elapsed.count() << " seconds" << std::endl;

    
  return 0;
}
