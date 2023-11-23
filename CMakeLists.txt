cmake_minimum_required(VERSION 3.10)

# Set some basic project attributes
project (SU2Chain
	VERSION 0.1
	DESCRIPTION "Su2 design optimization")

# This project will output an executable file
add_executable(${PROJECT_NAME} SU2Chain.cpp)

# Create a simple configuration header
configure_file(config.h.in config.h)

# Include the configuration header in the build
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")


find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIBRARIES})
