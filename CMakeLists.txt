cmake_minimum_required(VERSION 3.10)

# Set the project name
project(SM_KMP)

set(CMAKE_C_STANDARD 17)  # or 11, 17, etc.
set(CMAKE_C_STANDARD_REQUIRED True)

# Specify the path of the Gurobi solver
set(GUROBI_HOME /mnt/c/Users/Daniel/ubuntu/gurobi1201/linux64)
include_directories(${GUROBI_HOME}/include)
link_directories(${GUROBI_HOME}/lib)
set(GRB_LICENSE_FILE /mnt/c/Users/Daniel/ubuntu/gurobi1201/gurobi.lic)

#list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" ${CMAKE_MODULE_PATH})
#find_package(GUROBI REQUIRED)
#set(GUROBI_LIBRARY ${GUROBI_HOME}/lib/)
#set(GUROBI_INCLUDE_DIRS ${GUROBI_HOME}/include)
#message("-- GUROBI_LIBRARY points to " ${GUROBI_LIBRARY})
#message("-- GUROBI_INCLUDE_DIRS points to " ${GUROBI_INCLUDE_DIRS})

include_directories($ENV{GUROBI_HOME}/include)
link_directories($ENV{GUROBI_HOME}/lib)


# Specify the source directory
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

# Add the executable
add_executable(SMKMP 
    ${SOURCE_DIR}/macros.h
    ${SOURCE_DIR}/ts_config.h
    ${SOURCE_DIR}/item.h ${SOURCE_DIR}/item.c
    ${SOURCE_DIR}/ts_heur.h ${SOURCE_DIR}/ts_heur.c
    ${SOURCE_DIR}/data_structs.h ${SOURCE_DIR}/data_structs.c
    ${SOURCE_DIR}/cli_parse.h ${SOURCE_DIR}/cli_parse.c
    ${SOURCE_DIR}/instance.c ${SOURCE_DIR}/instance.h 
    ${SOURCE_DIR}/compactModel.cpp  ${SOURCE_DIR}/compactModel.h
    ${SOURCE_DIR}/labelSetting.cpp  ${SOURCE_DIR}/labelSetting.h
    ${SOURCE_DIR}/DualBound.cpp ${SOURCE_DIR}/ExactLinRelax.cpp
    ${SOURCE_DIR}/branch_and_pricing.cpp
    ${SOURCE_DIR}/binpackingHeuristic.cpp
    ${SOURCE_DIR}/main.cpp)
# Add the Gurobi library
target_link_libraries(SMKMP PRIVATE gurobi_c++ gurobi120 pthread)

# Set the output directories for object files and the executable
set_target_properties(SMKMP PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Ensure .o files are placed in the build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
# target_compile_options(knapsack PRIVATE -Wall -Wextra -Wpedantic -g -ggdb -O0)
target_compile_options(SMKMP PRIVATE -w -Wextra -Wpedantic -O3 -g -ggdb)
