cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME_STR go_multisolver)

project (${PROJECT_NAME_STR})

enable_language(C)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
find_package(OpenMP REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(3rd-party/global-optimization-test-problems)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/cmdline/
                    ${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/json)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

add_executable(${PROJECT_NAME_STR} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME_STR} gkls_functions grishagin_functions ${OpenMP_CXX_LIBRARIES})
set_target_properties(${PROJECT_NAME_STR} PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  target_compile_options(${PROJECT_NAME_STR} INTERFACE -Wall -march=native)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(/EHsc /openmp)
endif()
