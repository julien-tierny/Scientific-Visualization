# Copyright (C) Julien Tierny <julien.tierny@lip6.fr>

if(NOT WTFIT_PKG)

  # init project global variables
  set(PROJECT_DEP "" 
    CACHE INTERNAL "PROJECT_DEP")
  set(PROJECT_FLAGS "" 
    CACHE INTERNAL "PROJECT_FLAGS")
  set(PROJECT_SRC "" 
    CACHE INTERNAL "PROJECT_SRC")

  foreach(LIB ${LIB_LIST})
    set(${LIB} "" CACHE INTERNAL "${LIB}")
  endforeach(LIB)
  set(LIB_LIST "" CACHE INTERNAL "LIB_LIST")

  # specify the path to the core packages
  find_path(BASECODE_DIR baseCode.cmake
    PATHS
      ${WTFIT_DIR}/baseCode/
      baseCode/
      ../baseCode/
      ../../baseCode/
  )
  include(${BASECODE_DIR}/baseCode.cmake)

  # cmake helper functions
  function(wtfit_add_baseCode_package package)
    if(NOT ${package}_PATH)
    wtfit_find_package(${package})
      if(NOT ${package}_PATH)
        message(FATAL_ERROR 
          "Could not find the package '${package}.cmake'!") 
      endif(NOT ${package}_PATH)
    endif(NOT ${package}_PATH) 
  endfunction(wtfit_add_baseCode_package)

  function(wtfit_add_cflags flags)
    set(PROJECT_FLAGS ${PROJECT_FLAGS}
      ${flags}
      CACHE INTERNAL "PROJECT_FLAGS")
  endfunction(wtfit_add_cflags)

  function(wtfit_add_dep dep)
    set(PROJECT_DEP ${PROJECT_DEP} 
    ${dep} 
    CACHE INTERNAL "PROJECT_DEP")
  endfunction(wtfit_add_dep)

  function(wtfit_add_external_package header lib)
    option(with${lib} "Enable ${lib} support" true)

    if(with${lib})
      find_path(${header}_PATH ${header}
        PATHS
          ${ARGN}
      )
      if(NOT ${header}_PATH)
        message(STATUS
"wtfit -----------------------------------------------------------------------")
        message(STATUS "Header ${header} not found. External package disabled.")
        message(STATUS
"-----------------------------------------------------------------------------")
        set(with${lib} OFF)
      endif(NOT ${header}_PATH)

      if(${header}_PATH)
        include_directories(${${header}_PATH})
        find_library(${lib}_PATH 
          NAMES ${lib}
          PATHS
            ${ARGN}
        )
        if(NOT ${lib}_PATH)
          message(STATUS
"wtfit -----------------------------------------------------------------------")
          message(STATUS "Library ${lib} not found. External package disabled.")
          message(STATUS
"-----------------------------------------------------------------------------")
          set(with${lib} OFF)
        endif(NOT ${lib}_PATH)

        if(${lib}_PATH)
          wtfit_add_option(with${lib})
          wtfit_add_dep(${${lib}_PATH})
        endif(${lib}_PATH)
      endif(${header}_PATH)
    endif(with${lib})

    message(STATUS
"wtfit -----------------------------------------------------------------------")
    message(STATUS
    "${lib} support: ${with${lib}}   (-Dwith${lib}=)")
    message(STATUS
"-----------------------------------------------------------------------------")
  endfunction(wtfit_add_external_package)

  function(wtfit_add_option option)
    set(PROJECT_FLAGS "${PROJECT_FLAGS} -D${option}"
      CACHE INTERNAL "PROJECT_FLAGS")
  endfunction(wtfit_add_option)

  function(wtfit_add_optional_baseCode_package package)
    option(with${package} "Enable ${package} support" true)

    if(${with${package}})
      if(NOT ${package}_PATH)
        wtfit_find_package(${package})
        if(${package}_PATH)
          wtfit_add_option(with${package})
        elseif(NOT ${package}_PATH)
          set(with${package} OFF)
          message(STATUS
"wtfit -----------------------------------------------------------------------")
          message(STATUS
            "Package ${package} not found. Option disabled.")
          message(STATUS
"-----------------------------------------------------------------------------")
        endif(${package}_PATH)
      endif(NOT ${package}_PATH) 
    endif(${with${package}})

    message(STATUS
"wtfit -----------------------------------------------------------------------")
    message(STATUS
  "${package} support: ${with${package}}   (-Dwith${package}=)")
    message(STATUS
"-----------------------------------------------------------------------------")
  endfunction(wtfit_add_optional_baseCode_package)

  function(wtfit_add_optional_vtkWrapper_package package)
    if(NOT withVTK)
      message(FATAL_ERROR "Trying to compile a VTK wrapper with VTK disabled!")
    endif(NOT withVTK)
    wtfit_add_optional_baseCode_package(${package})
  endfunction(wtfit_add_optional_vtkWrapper_package)

  function(wtfit_add_source source)
    include_directories(${CMAKE_CURRENT_LIST_DIR})
    set(PROJECT_SRC ${PROJECT_SRC} 
      ${CMAKE_CURRENT_LIST_DIR}/${source}
      CACHE INTERNAL "PROJECT_SRC")
  endfunction(wtfit_add_source)

  function(wtfit_add_vtkWrapper_package package)
    if(NOT withVTK)
      message(FATAL_ERROR "Trying to compile a VTK wrapper with VTK disabled!")
    endif(NOT withVTK)
    wtfit_add_baseCode_package(${package})
  endfunction(wtfit_add_vtkWrapper_package)

  function(wtfit_find_package package)
    set(${package}_PATH ${package}_PATH-NOTFOUND 
      CACHE INTERNAL "${package}_PATH")
    find_path(${package}_PATH "${package}.cmake"
      PATHS
        ${CMAKE_BINARY_DIR}/../../../../sandbox/baseCode/${package}/
        ${CMAKE_BINARY_DIR}/../../../../sandbox/vtkWrappers/${package}/
        ${CMAKE_BINARY_DIR}/../../../sandbox/baseCode/${package}/
        ${CMAKE_BINARY_DIR}/../../../sandbox/vtkWrappers/${package}/
        ${CMAKE_BINARY_DIR}/../../sandbox/baseCode/${package}/
        ${CMAKE_BINARY_DIR}/../../sandbox/vtkWrappers/${package}/
        ${CMAKE_BINARY_DIR}/../sandbox/baseCode/${package}/
        ${CMAKE_BINARY_DIR}/../sandbox/vtkWrappers/${package}/
        ${CMAKE_BINARY_DIR}/sandbox/baseCode/${package}/
        ${CMAKE_BINARY_DIR}/sandbox/vtkWrappers/${package}/
        ${WTFIT_DIR}/baseCode/${package}/
        ${WTFIT_DIR}/vtkWrappers/${package}/
      NO_DEFAULT_PATH)
    if(${package}_PATH)
      # the package may just be a template (no call to add_source)
      include_directories("${${package}_PATH}")
      include("${${package}_PATH}/${package}.cmake")
    endif(${package}_PATH)
  endfunction(wtfit_find_package)

  function(wtfit_wrapup_binary project)
    if(APPLE)
      add_executable(${project} MACOSX_BUNDLE ${PROJECT_SRC})
    else()
      add_executable(${project} ${PROJECT_SRC})
    endif(APPLE)
    wtfit_wrapup_flags(${project})
  endfunction(wtfit_wrapup_binary)

  function(wtfit_wrapup_flags project)
    # set up compilation flags pulled out by the called modules 
    set_target_properties(${project}
      PROPERTIES
      COMPILE_FLAGS 
      "${PROJECT_FLAGS}")
  
    # specify the required libraries for linkage, as pulled out by the called 
    # modules (${PROJECT_DEP})
    target_link_libraries(${project} "${PROJECT_DEP}")

    message(STATUS
"wtfit -----------------------------------------------------------------------")
    message(STATUS "wtfit project: ${project}")
    get_property(PROJECT_INCLUDES DIRECTORY 
      ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
    message(STATUS "wtfit includes: ${PROJECT_INCLUDES}")
    message(STATUS "wtfit deps: ${PROJECT_DEP}")
    message(STATUS "wtfit flags: ${PROJECT_FLAGS}")
    message(STATUS "wtfit src: ${PROJECT_SRC}")
    message(STATUS
"-----------------------------------------------------------------------------")
  endfunction(wtfit_wrapup_flags)

  function(wtfit_wrapup_library library source)

    # link only once!
    foreach(LIB ${LIB_LIST})
      if("${LIB}" MATCHES " ${source}")
        return()
      endif("${LIB}" MATCHES " ${source}")
    endforeach(LIB)

    foreach(sourceFile "${source}")
      add_library(${library} OBJECT
        ${CMAKE_CURRENT_LIST_DIR}/${sourceFile})
    endforeach(sourceFile)


    # use the project flags
    set_target_properties(${library}
      PROPERTIES 
      COMPILE_FLAGS 
      "${PROJECT_FLAGS}")

    # add the library to the calling project
    add_library(${PROJECT_NAME} SHARED $<TARGET_OBJECTS:${library}>)

    set(LIB_LIST "${LIB_LIST} ${source}" CACHE INTERNAL "LIB_LIST")

  endfunction(wtfit_wrapup_library)

  function(wtfit_wrapup_paraview_plugin
    plugin_name
    plugin_version)

    if(withParaView)
      
      add_paraview_plugin(
        # name
        ${plugin_name}
        # version
        "${plugin_version}"                         
        SERVER_MANAGER_XML "${plugin_name}.xml"
        SERVER_MANAGER_SOURCES "${PROJECT_SRC}"
      )
    endif(withParaView)

    wtfit_wrapup_flags(${plugin_name})
  endfunction(wtfit_wrapup_paraview_plugin)

  # message about the location of the common code base
  message(STATUS
"wtfit -----------------------------------------------------------------------")
  message(STATUS "wtfit path: ${WTFIT_DIR}")
  message(STATUS
"-----------------------------------------------------------------------------")

  option(withVTK "Enable VTK support" true)

  option(withParaView "Enable ParaView support" false)

  option(withVTKGUI "Enable VTK GUI support" false)
  
  # add the common package
  wtfit_add_baseCode_package(common)
  
  if(NOT withParaView)
    if(withVTK)
      if(NOT VTK_USE_FILE)
        find_package(VTK REQUIRED)
        include(${VTK_USE_FILE})  
        set(PROJECT_DEP "${PROJECT_DEP} ${VTK_LIBRARIES}" 
          CACHE INTERNAL "PROJECT_DEP")
        set(PROJECT_FLAGS "${PROJECT_FLAGS} -DwithVTK"
          CACHE INTERNAL "PROJECT_FLAGS")
        if(withVTKGUI)
          wtfit_add_vtkWrapper_package(vtkTextureMapFromField)
          wtfit_add_vtkWrapper_package(vtkWRLExporter)
        endif(withVTKGUI)
      endif(NOT VTK_USE_FILE)
    endif(withVTK)
  endif(NOT withParaView)

  if(withParaView)
    find_package(ParaView REQUIRED)
    include(${PARAVIEW_USE_FILE})
  endif(withParaView)

  # messages on dependencies
  message(STATUS
"wtfit -----------------------------------------------------------------------")
  message(STATUS "baseCode path: ${BASECODE_DIR}")

  # optional dependency on VTK
  message(STATUS "VTK support: ${withVTK}   (-DwithVTK=)")
  message(STATUS "VTK GUI support: ${withVTKGUI}   (-DwithVTKGUI=)")
  message(STATUS "ParaView support: ${withParaView}   (-DwithParaView=)")
  message(STATUS
"-----------------------------------------------------------------------------")

  set(WTFIT_PKG true)

endif(NOT WTFIT_PKG)
