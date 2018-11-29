cmake_minimum_required(VERSION 2.8)

####
# project settings
####

project(RayTrace)

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
	raytrace.cpp
)


####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)

include_directories(${GLAMER_INCLUDE_DIRS})


####
# build
####

add_executable(RayTrace ${HEADERS} ${SOURCES})

target_link_libraries(RayTrace ${GLAMER_LIBRARIES})


