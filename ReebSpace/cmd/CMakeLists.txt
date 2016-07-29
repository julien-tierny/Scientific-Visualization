cmake_minimum_required(VERSION 2.4)

# name of the project (dummy)
project(reebSpaceCmd)

# path to the common code base
find_path(WTFIT_DIR wtfit.cmake
  PATHS
    wtfit/
    ../wtfit/
    ../../wtfit/
    ../../../wtfit/
)
include(${WTFIT_DIR}/wtfit.cmake)

# include the different packages required for the application
wtfit_add_vtkWrapper_package(vtkReebSpace)

# add the local cpp files 
wtfit_add_source(Editor.cpp)
wtfit_add_source(main.cpp)

# wrap-up the build
wtfit_wrapup_binary(reebSpaceCmd)
