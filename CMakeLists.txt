cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME_STR go_multisolver)

project(${PROJECT_NAME_STR})

enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
find_package(OpenMP REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(3rd-party/global-optimization-test-problems)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/cmdline/
                    ${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/json)

add_library(multisolver_core ${CMAKE_CURRENT_SOURCE_DIR}/src/evolvent.cpp)
target_include_directories(multisolver_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(solve_sample ${CMAKE_CURRENT_SOURCE_DIR}/src/solve_sample.cpp)
target_link_libraries(solve_sample multisolver_core gkls_functions grishagin_functions ${OpenMP_CXX_LIBRARIES})
set_target_properties(solve_sample PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})

add_executable(solve_constrained_sample ${CMAKE_CURRENT_SOURCE_DIR}/src/solve_constrained_sample.cpp)
target_link_libraries(solve_constrained_sample multisolver_core ${OpenMP_CXX_LIBRARIES})
set_target_properties(solve_constrained_sample PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(solve_sample INTERFACE -Wall -march=native)
  target_compile_options(solve_constrained_sample INTERFACE -Wall -march=native)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(/EHsc /openmp)
endif()
