cmake_minimum_required(VERSION 3.12)
project(KokkosSundialODESolver CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Kokkos and ADIOS2 locations
set(Kokkos_DIR "$ENV{HOME}/kokkos-install/lib/cmake/Kokkos")
set(CMAKE_BUILD_TYPE Debug)

# Find packages (Sundials manually)
find_package(Boost REQUIRED)
find_package(Kokkos REQUIRED)
find_package(adios2 REQUIRED)

# Include directories
include_directories(
    ${Boost_INCLUDE_DIRS}
    $ENV{HOME}/kokkos-install/include
    $ENV{HOME}/sundials-install/include
)

# Link directories
link_directories(
    $ENV{HOME}/kokkos-install/lib
    $ENV{HOME}/sundials-install/lib
)

# Solver library (uses Sundials)
add_library(SolverLibrary STATIC
    SUNDIAL.cpp
    boost.cpp
    rungekutta.cpp
)

# Define for static SUNDIALS linking if needed
target_compile_definitions(SolverLibrary PRIVATE SUNDIALS_STATIC)

# Link Sundials to SolverLibrary
target_link_libraries(SolverLibrary
    sundials_cvode
    sundials_nvecserial
    sundials_sunmatrixdense
    sundials_sunlinsoldense
    sundials_core
    m
    dl
)

# Main executable
add_executable(Electromechanical_system_kokkos
    Electromechanical_system_kokkos.cpp
)
add_executable(Electromechanical_system
    Electromechanical_system.cpp
)

# Link everything together
# target_link_libraries(Electromechanical_system_kokkos
#     SolverLibrary
#     Kokkos::kokkos
#     adios2::adios2
# )
target_link_libraries(Electromechanical_system
    SolverLibrary
    Kokkos::kokkos
    adios2::adios2
)
