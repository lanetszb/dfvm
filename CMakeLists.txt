cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_BUILD_TYPE Release)

set(POJECT_NAME dfvm)
project(${POJECT_NAME})
add_subdirectory(sgrid)
add_subdirectory(diffusion)

