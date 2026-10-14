cmake_minimum_required(VERSION 3.16)

project(StochasticDiffusion LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Ceres REQUIRED)
find_package(Boost CONFIG REQUIRED)

add_executable(stochastic_diffusion
    src/main.cpp
    src/io_handler.cpp
    src/model_data.cpp
    src/model_functions.cpp
    src/visualization.cpp
    src/regression.cpp
    src/estimations.cpp
    src/equations.cpp
    src/optimization.cpp
    src/evolutive.cpp
)

target_include_directories(stochastic_diffusion PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(stochastic_diffusion PRIVATE
    Eigen3::Eigen
    Ceres::ceres
    Boost::boost
)

target_compile_options(stochastic_diffusion PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
)

enable_testing()

add_test(
    NAME smoke
    COMMAND stochastic_diffusion
)

set_tests_properties(
    smoke
    PROPERTIES
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        TIMEOUT 10
)
