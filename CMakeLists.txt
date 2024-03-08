cmake_minimum_required(VERSION 3.0...3.8)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#set(CMAKE_C_COMPILER "/usr/local/Cellar/gcc/13.2.0/bin/gcc-13")
#set(CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc/13.2.0/bin/c++-13")

####
# project settings
####

project(example)

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


####
# headers
####

set(HEADERS
)


####
# sources
####

set(SOURCES
	example.cpp
)


####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)

include_directories(${GLAMER_INCLUDE_DIRS})


####
# build
####

add_executable(example ${HEADERS} ${SOURCES})

target_link_libraries(example ${GLAMER_LIBRARIES})


