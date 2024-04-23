cmake_minimum_required(VERSION 3.5)

project(bridge)

set(CMAKE_CXX_FLAGS "${CXXFLAGS}")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(GUROBI REQUIRED)

if(NOT GUROBI_FOUND)
    message(SEND_ERROR "couldn't find Gurobi")
else()
    message("found Gurobi")
    cmake_print_variables(GUROBI_INCLUDE_DIRS)
    include_directories(${GUROBI_INCLUDE_DIRS})
endif()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(polyscope)

add_library(triangle triangle/triangle.c)

add_executable(${CMAKE_PROJECT_NAME} bridge.cpp)

target_link_libraries(${PROJECT_NAME} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

target_link_libraries (${PROJECT_NAME} Eigen3::Eigen)

target_link_libraries(${PROJECT_NAME} polyscope triangle)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)


