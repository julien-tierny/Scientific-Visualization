if(NOT SURFACE_SCALAR_FIELD_PKG)
  
  find_path(MESH_PATH mesh.cmake
    PATHS
    baseCode/mesh/)

  if(MESH_PATH)
    # depends on the mesh module
    include(baseCode/mesh/mesh.cmake) 
  elseif(NOT MESH_PATH)
    message (FATAL_ERROR "Could not find the mesh package! (mesh.cmake)")
  endif(MESH_PATH)

  include_directories("baseCode/common/" "baseCode/surfaceScalarField")

  set(PROJECT_SRC ${PROJECT_SRC}  
    baseCode/surfaceScalarField/SurfaceScalarField.cpp)

  set(SURFACE_SCALAR_FIELD_PKG true)
endif(NOT SURFACE_SCALAR_FIELD_PKG)
