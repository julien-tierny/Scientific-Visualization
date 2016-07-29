wtfit_add_baseCode_package(zeroSkeleton)

wtfit_add_vtkWrapper_package(vtkOneSkeleton)

# This class is not meant to be used by a paraview plugin, use library.
wtfit_wrapup_library(libvtkZeroSkeleton "vtkZeroSkeleton.cpp")