cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE "Release")

project(fiberSurface)

option(withBVH "Enable BVH support" false)

if(withBVH)
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DwithBVH")
  set(CMAKE_CXX_FLAGS 
    "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DwithBVH")
endif(withBVH)

option(withOpenMP "Enable OpenMP support" true)
if(withOpenMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS 
      "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DwithOpenMP")
    set(CMAKE_CXX_FLAGS 
      "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DwithOpenMP")
  endif(OPENMP_FOUND)
endif(withOpenMP)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories("bvh")
include_directories("fiberSurface")
include_directories("octree")
include_directories("wtfit")

# build bvh library
add_library(bvh SHARED bvh/bvh.c)

set_target_properties(bvh
  PROPERTIES
  COMPILE_FLAGS "-march=native -O3 -std=c99 ${CMAKE_C_FLAGS} -mavx -DNDEBUG")

# build fiber surface library
add_library(FiberSurface SHARED 
  fiberSurface/vtkFiberSurfaceFilter.cpp
  fiberSurface/vtkCleanUnstructuredGrid.cxx)
set_target_properties(FiberSurface
  PROPERTIES
  COMPILE_FLAGS "-march=native -O3 -std=c++11 ${CMAKE_CXX_FLAGS}")
  
# build wtfit
add_library(wtfit SHARED wtfit/Debug.cpp)
set_target_properties(wtfit
  PROPERTIES
  COMPILE_FLAGS "-march=native -O3 -std=c++11 ${CMAKE_CXX_FLAGS}")

# build binary
add_executable(fiberSurface
    main.cpp)

set_target_properties(fiberSurface
  PROPERTIES
  COMPILE_FLAGS "-std=c++11 -march=native -O3 ${CMAKE_CXX_FLAGS}")

target_link_libraries(fiberSurface 
  "${VTK_LIBRARIES}" bvh wtfit FiberSurface)
