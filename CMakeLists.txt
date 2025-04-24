cmake_minimum_required(VERSION 3.12)
project(KokkosSundialODESolver CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Kokkos and ADIOS2 locations
set(Kokkos_DIR "$ENV{HOME}/kokkos-install/lib/cmake/Kokkos")
set(CMAKE_BUILD_TYPE Debug)

# Find packages (Sundials manually)
find_package(Boost REQUIRED)
find_package(Kokkos REQUIRED)
find_package(ADIOS2 REQUIRED)

# Include directories
include_directories(
    ${Boost_INCLUDE_DIRS}
    $ENV{HOME}/kokkos-install/include
    $ENV{HOME}/sundials-install/include
)
add_subdirectory(externals/Catch2)

# Link directories
link_directories(
    $ENV{HOME}/kokkos-install/lib
    $ENV{HOME}/sundials-install/lib
)

# Solver library (uses Sundials)
add_library(SolverLibrary STATIC
    Solvers/SUNDIAL.cpp
    Solvers/boost.cpp
    Solvers/rungekutta.cpp
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
add_executable(Electromechanical_system
    Electromechanical_system.cpp
)
# Catch2 tests
add_executable(test_sundials Tests/test_sundial.cpp Solvers/SolverLibrary.hpp)
add_executable(test_rk4 Tests/test_rk4.cpp Solvers/SolverLibrary.hpp)
add_executable(test_adios Tests/test_adios.cpp)

# Add include directories
target_include_directories(test_sundials PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(test_adios PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(test_rk4 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Link required Sundials libraries
target_link_libraries(test_sundials
  PUBLIC 
    SolverLibrary
  PRIVATE
    Catch2::Catch2WithMain
    sundials_cvode
    sundials_nvecserial
    sundials_sunlinsoldense
    sundials_sunmatrixdense
    dl
)
target_link_libraries(Electromechanical_system
     SolverLibrary
     Kokkos::kokkos
     adios2::adios2
)
target_link_libraries(test_adios PRIVATE adios2::core adios2::cxx11 Catch2::Catch2WithMain)
target_link_libraries(test_rk4 PUBLIC SolverLibrary PRIVATE Catch2::Catch2WithMain)

# Register with CTest
include(CTest)
include(Catch)
catch_discover_tests(test_sundials)
catch_discover_tests(test_adios)
catch_discover_tests(test_rk4)


# add_executable(Electromechanical_system_kokkos
#     Electromechanical_system_kokkos.cpp
# )
# Link everything together
# target_link_libraries(Electromechanical_system_kokkos
#     SolverLibrary
#     Kokkos::kokkos
#     adios2::adios2
# )




