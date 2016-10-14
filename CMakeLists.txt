cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME_STR go_miltisolver)

project (${PROJECT_NAME_STR})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  add_definitions(-std=c++11 -Wall -fopenmp)
elseif (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  add_definitions(/EHsc /openmp)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3d-party)

file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/GKLS/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/GKLS/*.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Grishagin/*.hpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/Grishagin/*.cpp)

add_executable(${PROJECT_NAME_STR} ${SRC_FILES})
set_target_properties(${PROJECT_NAME_STR} PROPERTIES LINK_FLAGS "-fopenmp")
