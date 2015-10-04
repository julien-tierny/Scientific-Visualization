if(NOT MESH_PKG)
  
  include_directories("baseCode/common/" "baseCode/mesh")

  set(PROJECT_SRC ${PROJECT_SRC} 
    baseCode/mesh/TriangleMesh.cpp)

  set(PROJECT_FLAGS ${PROJECT_FLAGS} "-O3")

  set(MESH_PKG true)
endif(NOT MESH_PKG)
