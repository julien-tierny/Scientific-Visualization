cmake_minimum_required(VERSION 2.4)
PROJECT(generalSimplification)


#add the different modules required by the application
include(baseCode/surfaceScalarField/surfaceScalarField.cmake)

add_executable(generalSimplification ${PROJECT_SRC} 
  Editor.cpp
  main.cpp)

set_target_properties(generalSimplification
  PROPERTIES
  COMPILE_FLAGS
  "${PROJECT_FLAGS}")

target_link_libraries(generalSimplification ${PROJECT_DEP})
