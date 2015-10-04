DATA DISCLAIMER
---------------
The data sets contained in this archive are the property of their copyright
owners. They are made available in this archive for review purpose only. They 
should not be redistributed publicly without their owner's permission.

CODE DISCLAIMER
---------------
This code is provided "as is" without warranty of any kind. It is made available
for review purpose only.




INSTRUCTIONS
------------
The following instructions are for Unix environments (MacOS, Linux, etc.). For
Windows users, please follow the usual procedure to import a CMake project into
VisualStudio (please be careful though, important directory paths may need to 
be changed for the program to run properly).


IMPORTANT NOTE
--------------
Due to submission storage limitations (100MB), only a toy example data-set is 
provided in this archive.



1) Software requirements:
-------------------------
This code requires the VTK development libraries to be installed on the system 
(version 6.1 or higher). Please see http://www.vtk.org/ for installation 
instructions.
A recent version of CMake (http://www.cmake.org) should be installed (version 
2.8 or higher).


2) Compiling the code:
----------------------
After decompressing the archive, enter the following directory:

fiberSurface/src/build

From there, enter the following commands (omit the "$" character):

$ cmake ../
$ make


IMPORTANT NOTE
--------------
Our BVH implementation uses AVX features. If your CPU was purchased BEFORE 
2011, it is very likely it does not support them.
Thus, by default, our BVH implementation is disabled at build time. With the 
default build, our code will default back to the octree algorithm if you select 
the BVH algorithm, and it will display the following error message:

[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter] Not compiled with BVH support!
[vtkFiberSurfaceFilter] If your CPU supports AVX features, re-run cmake with the 
option '-DwithBVH=ON'
[vtkFiberSurfaceFilter] Defaulting back to the octree implementation...
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]
[vtkFiberSurfaceFilter]

If your CPU supports AVX features, you can enable our BVH support with the 
following command lines (from the directory "fiberSurface/src/build"):

$ cmake ../ -DwithBVH=ON
$ make

Note that if you enable the BVH support while your CPU does not support AVX 
features, our implementation is likely to crash if the BVH algorithm is 
selected.


3) Running the code:
--------------------
The directory "fiberSurface/data/" contains a bivariate field defined on a 
tetrahedral mesh that has been obtained by a 5-subdivision of a DOWNSAMPLED 
version of our EthaneDiol data set.
Therefore, due to submission storage limitations (100MB), you will not be able 
to reproduce exactly the result illustrated in the paper.

After building the code, enter the following directory:

fiberSurface/src/build

And run the following command lines:

./fiberSurface \
  -i ../../data/ethaneDiol_downSampled.vtu \
  -iF0 1 -iF1 0 \
  -p ../../data/teaser_polygon.vtu

./fiberSurface \
  -i ../../data/ethaneDiol_downSampled.vtu \
  -iF0 1 -iF1 0 \
  -p ../../data/black_polygon.vtu

./fiberSurface \
  -i ../../data/ethaneDiol_downSampled.vtu \
  -iF0 1 -iF1 0 \
  -p ../../data/blue_polygon.vtu

./fiberSurface \
  -i ../../data/ethaneDiol_downSampled.vtu \
  -iF0 1 -iF1 0 \
  -p ../../data/green_polygon.vtu

./fiberSurface \
  -i ../../data/ethaneDiol_downSampled.vtu \
  -iF0 1 -iF1 0 \
  -p ../../data/red_polygon.vtu


At each execution of the program "fiberSurface", the output fiber surface will 
be save in a file named "output.vtp". 
Thus, you may want to rename this file in between each execution.

For each execution of the program, the console output should be of that effect:
[fiberSurface] Reading input data-set...
[fiberSurface]   112995 tet(s) read.
[fiberSurface] Reading input polygon...
[fiberSurface]   5 segment(s) read.
[fiberSurface] Data fields: 'log(Rho)' and 'log(s)'
[fiberSurface] Polygon fields: 'u' and 'v'
[fiberSurface] Implementation: 0
[fiberSurface] Algorithm: 2
[fiberSurface] Thread number: 1
[fiberSurface] Threading strategy: 1
[fiberSurface] Octree minimum range area ratio: 0
[fiberSurface] BVH minimum tet number: 1
[fiberSurface] Automatic threading strategy balance: 1
[fiberSurface] Manifold output: 0
[fiberSurface] Fiber texture coordinates: 0
[fiberSurface] Polygon edge segmentation: 0
[fiberSurface] Let's go!
[fiberSurface]
[vtkFiberSurfaceFilter] Using 1 thread(s)...
[vtkFiberSurfaceFilter] Processing polygon edge #0
[vtkFiberSurfaceFilter] Processing polygon edge #1
[vtkFiberSurfaceFilter] Processing polygon edge #2
[vtkFiberSurfaceFilter] Processing polygon edge #3
[vtkFiberSurfaceFilter] Processing polygon edge #4
[vtkFiberSurfaceFilter] Reduce step done in 0.000658035 s.
[vtkFiberSurfaceFilter] Fiber surface computed in 0.017894 s. TOTAL (47724 #v, 
15908 #f)
[fiberSurface] Saving output to output.vtp...
[fiberSurface]   done!

The output fiber surface is saved in a VTK surface file format that can be read 
directly from ParaView.
At this point, if you load all the generated surfaces, you should obtain a 
similar visualization to the one illustrated in the following file:
"fiberSurface/data/screenShot.png"


4) Program options and arguments:
---------------------------------
From the directory "fiberSurface/src/build", running the following command line 
will display the program options and arguments:
./fiberSurface -h

Of particular interest are the following optional arguments:
-t <Thread number (default: 1)>
-I <Implementation: [0, Regular] [1, Octree] [2, BVH] (default: 0)>


4) Using other data-sets:
-------------------------
The input bivariate tetrahedral mesh must be saved as a vtkUnstructuredGrid 
(VTU file format). Note that only pure tet-meshes are supported. Make sure your 
geometry is tet-only.

Each component of the bivariate field should be attached as an individual 
scalar field to the vertices of the geometry.
To select the scalar fields to consider as components of your bivariate field, 
use the following arguments:
-iF0 <Scalar field Id for first component (data-set, default: 0)>
-iF1 <Scalar field Id for second component (data-set, default: 1)>

The input fiber surface control polygon (FSCP) must be saved as a 
vtkUnstructuredGrid (VTU file format) as a list of edges.
The range coordinates of each FSCP vertex should be expressed through 
individual scalar fields attached to the vertices of the FSCP geometry.
To select the scalar fields to consider as range coordinates of your FSCP, use 
the following arguments:
-pF0 <Scalar field Id for first component (polygon, default: 0)>
-pF1 <Scalar field Id for second component (polygon, default: 1)>