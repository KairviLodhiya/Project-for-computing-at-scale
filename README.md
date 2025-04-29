# Project for Computing at Scale

---  

*This project implements an electromechanical system solver using C++, Kokkos for performance portability, Sundials (CVODE) for ODE integration, and ADIOS2 for data I/O. Unit tests are implemented using Catch2.*

---

### Need: 
#### This code is used in solving system of equation for a non-linear electro-mechanical coupled oscillator. This oscillator mimics a part of an low powered sensor. This sensor is understood to improve the energy transfer from electrical part to the mechanical part of the system. The system of equations are as follows:
---
                                                    my" + ky - 1/2 q^2/eA = kg
                                                    Lq" + q/C0 + qy/eA = 0
---
#### The Non-dimensional version of this is:
---
                                                    y" + y + 1/4 C1q^2 = 1/2
                                                    q" + C2q + 2C3qy = 0
---                                                    
#### Here C1, C2 and C3 are constant parameters which define the system's internal relations. The goal is the toggle these parameters and the boundary conditions [y0 ,q0 ,y0',q0'] to see the nature of the system. These parameters have internal relation - C2+C3 = 1 to equate natural frequencies of both the oscillators. 

---

### Inputs and Outputs of the software

- **Inputs**: The input would be two equations and boundary conditions. It would also ask for number of time steps that this system needs to be running.
  This can be done by updating the equations and passing the initial conditions in the function called cv_derivatives() 

- **Outputs**: The output of the system would be a graphs for y(t), q(t), y'(t) and q'(t) and saves a file with numerical solution of the equation for the number of steps given.

---

### Steps to build and run the Project

- [Kokkos](https://github.com/kokkos/kokkos) — for performance portability
- [Sundials (CVODE)](https://computing.llnl.gov/projects/sundials) — for ODE integration
- [ADIOS2](https://adios2.readthedocs.io/) — for high-performance I/O
- [Catch2](https://github.com/catchorg/Catch2) — for unit testing
- [Python](https://www.python.org/) (for post-processing)
  - Packages: `numpy`, `matplotlib`, `mpi4py`
- MPI (for parallel visualization support)


- **Pre-requisites**

  - Linux or WSL (Ubuntu 20.04+ recommended)
  - GCC 9+, Clang 10+, or another C++17-capable compiler
  - CMake ≥ 3.16
  - Python 3.x (for post-processing)
  - Git with submodule support

#### Steps:
  1. - **Clone the Repo**

    git clone --recurse-submodules https://github.com/KairviLodhiya/Project-for-computing-at-scale.git
    cd Project-for-computing-at-scale


  2. - **Build and Install Kokkos**

    git clone https://github.com/kokkos/kokkos.git
    mkdir -p kokkos/build
    cd kokkos/build
    cmake .. \
      -DCMAKE_CXX_COMPILER=$HOME/openmpi/bin/mpic++ \
      -DKokkos_ENABLE_OPENMP=ON \
      -DKokkos_ENABLE_SERIAL=ON \
      -DCMAKE_INSTALL_PREFIX=$HOME/kokkos-install
    make -j$(nproc)
    make install
    cd ../..


  3. - **Build and install Sundials from Submodule**

    cd sundials
    mkdir build && cd build
    cmake .. \
      -DCMAKE_INSTALL_PREFIX=$HOME/sundials-install \
      -DBUILD_SHARED_LIBS=ON \
      -DSUNDIALS_ENABLE_CXX=ON \
      -DSUNDIALS_ENABLE_EXAMPLES=OFF \
      -DSUNDIALS_ENABLE_TESTS=OFF
    make -j$(nproc)
    make install
    cd ../../..


  4. - **Build and Install ADIOS2**

    cd ADIOS2
    mkdir build && cd build
    cmake .. \
      -DADIOS2_USE_MPI=OFF \
      -DADIOS2_USE_Python=ON \
      -DADIOS2_BUILD_BINDINGS=ON \
      -DCMAKE_INSTALL_PREFIX=$HOME/adios2-install
    make -j$(nproc)
    make install
    cd ../../..

  5. - **Build the project**

    mkdir build
    cd build
    cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=$HOME/openmpi/bin/mpic++ \
      -DCMAKE_PREFIX_PATH="$HOME/kokkos-install;$HOME/adios2-install;$HOME/sundials-install" \
      -DADIOS2_DIR=$HOME/adios2-install/lib/cmake/adios2
    make -j$(nproc)

  
  6. - **Run the executable**

    ./Electromechanical_system

  7. - **Run the unit tests**

    ctest --output-on-failure

- **Post-Processing** 

    Plot the results from the ADIOS file saved. 
    export PYTHONPATH=$HOME/adios2-install/lib/python3.10/site-packages:$PYTHONPATH
    python3 plot_adios.ipynb



