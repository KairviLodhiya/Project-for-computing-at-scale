#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>
#include <adios2.h>
#include <vector>

typedef double real;
constexpr int STATE_SIZE = 5; // time, y, q, dy, dq

// Implementing GPU based RK4
struct ElectromechanicalSystem {
  double C1, C2, C3;

  KOKKOS_INLINE_FUNCTION
  ElectromechanicalSystem(double C1_in = 0.2, double C2_in = 0.7)
      : C1(C1_in), C2(C2_in), C3(1.0 - C2_in) {}

  // y = [y, y_dot, q, q_dot]
  KOKKOS_INLINE_FUNCTION
  void operator()(const real* y, real* dydt) const {
    dydt[0] = y[2]; 
    dydt[1] = y[3];
    dydt[2] = 0.5 - y[0] - (C1 / 4.0) * y[1] * y[1]; // dq/dt = q'
    dydt[3] = -C2 * y[1] - 2.0 * C3 * y[1] * y[0];
  }
};

// RK4 step
KOKKOS_INLINE_FUNCTION
void rk4_step(const ElectromechanicalSystem& sys, real* y, real dt) {
  real dydt[4], k1[4], k2[4], k3[4], k4[4], y_temp[4];

  sys(y, k1);
  for (int i = 0; i < 4; ++i) y_temp[i] = y[i] + 0.5 * dt * k1[i];
  sys(y_temp, k2);
  for (int i = 0; i < 4; ++i) y_temp[i] = y[i] + 0.5 * dt * k2[i];
  sys(y_temp, k3);
  for (int i = 0; i < 4; ++i) y_temp[i] = y[i] + dt * k3[i];
  sys(y_temp, k4);
  for (int i = 0; i < 4; ++i)
    y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    const int N = 10;
    const real t_start = 0.0;
    const real t_end = 5000.0;   // for test
    const real dt = 0.001;
    const int steps = static_cast<int>((t_end - t_start) / dt) + 1;

    // 3D View: [simulation][step][state]
    Kokkos::View<real***> results("results", N, steps, STATE_SIZE);

    // Run all simulations in parallel on GPU
    Kokkos::parallel_for("RK4Sim", Kokkos::RangePolicy<>(0, N), KOKKOS_LAMBDA(int i) {
      ElectromechanicalSystem sys;
      real y[4];
      y[0] = i * 0.1;   // y
      y[1] = 1.0;       // q
      y[2] = 0.0;       // dy
      y[3] = 0.0;       // dq

      real t = t_start;
      for (int step = 0; step < steps; ++step) {
        results(i, step, 0) = t;
        results(i, step, 1) = y[0];  // y
        results(i, step, 2) = y[1];  // q
        results(i, step, 3) = y[2];  // dy
        results(i, step, 4) = y[3];  // dq
        rk4_step(sys, y, dt);
        t += dt;
      }
    });

    // Copy results to host for ADIOS2 output
    auto results_host = Kokkos::create_mirror_view(results);
    Kokkos::deep_copy(results_host, results);

    // ADIOS2 write
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("io");
    adios2::Engine writer = io.Open("gpu_parallel_rk4.bp", adios2::Mode::Write);

    for (int i = 0; i < N; ++i) {
      std::string group = "y" + std::to_string(static_cast<int>(i * 10)) + "_q100";

      std::vector<real> time(steps), y(steps), q(steps), dy(steps), dq(steps);
      for (int step = 0; step < steps; ++step) {
        time[step] = results_host(i, step, 0);
        y[step]    = results_host(i, step, 1);
        q[step]    = results_host(i, step, 2);
        dy[step]   = results_host(i, step, 3);
        dq[step]   = results_host(i, step, 4);
      }

      auto v_time = io.DefineVariable<real>(group + "/time", {}, {}, {steps});
      auto v_y    = io.DefineVariable<real>(group + "/y", {}, {}, {steps});
      auto v_q    = io.DefineVariable<real>(group + "/q", {}, {}, {steps});
      auto v_dy   = io.DefineVariable<real>(group + "/dy", {}, {}, {steps});
      auto v_dq   = io.DefineVariable<real>(group + "/dq", {}, {}, {steps});

      writer.BeginStep();
      writer.Put(v_time, time.data());
      writer.Put(v_y, y.data());
      writer.Put(v_q, q.data());
      writer.Put(v_dy, dy.data());
      writer.Put(v_dq, dq.data());
      writer.EndStep();
    }

    writer.Close();
    std::cout << "[DONE] All simulations completed and saved." << std::endl;
  }
  Kokkos::finalize();
  return 0;
}
