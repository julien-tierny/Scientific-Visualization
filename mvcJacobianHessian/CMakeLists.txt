cmake_minimum_required(VERSION 2.4)
PROJECT(mvcTest)

SET(CMAKE_BUILD_TYPE "Release")

include_directories("include" "shared")

# specify an executable 
add_executable(mvcTest 
  main.cpp)

