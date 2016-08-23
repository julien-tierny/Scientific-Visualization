/// \ingroup vtkWrappers
/// \class vtkTriangulation
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date January 2016.
/// 
/// \brief vtkTriangulation is a class that provides time and memory efficient 
/// traversal methods on triangulations of piecewise linear manifolds. It 
/// provides the following features:
///   -# Given a vertex, it provides: the list of edges that are connected to 
/// it, the list of its neighbors, its link, its star, etc.
///   -# Given an edge, it provides: its vertices, its star, etc.
///   -# Given a triangle, its provides: its vertices, its edges, etc.
///   -# Given a tetrahedron, its provides: its vertices, its edges, its 
/// neighbor tetrahedra, etc.
///   -# Given a triangulation, it provides: its list of vertices, edges, 
/// triangles and tetrahedra.
///
/// vtkTriangulation implements faster accesses than, more general-purpose, 
/// competing VTK data structures such as the vtkUnstructuredGrid class.
/// 
/// vtkTriangulation supports both explicit and implicit triangulations: 
///   -# Explicit triangulations with vtkUnstructuredGrid or vtkPolyData 
/// objects: Given a vtkUnstructuredGrid or a vtkPolyData representing a valid 
/// triangulation, vtkTriangulation provides time efficient accesses (requiring 
/// adequate pre-processing, see the Triangulation class documentation).
///   -# Implicit triangulations with vtkImageData objects: Given a vtkImageData
/// representing a regular grid, vtkTriangulation will perform an implicit 
/// triangulation of the grid, enabling both time and memory efficient 
/// traversals of triangulations of regular grids.
///
/// Apart from pre-processes, vtkTriangulation requires no memory overhead in
/// addition to the input vtkUnstructuredGrid, vtkPolyData or vtkImageData 
/// objects (the actual data is passed through pointers).
///
/// \note
/// Only pre-process the information you need! See the 
/// wtfit::Triangulation class documentation.
/// \sa wtfit::Triangulation

#ifndef _VTK_TRIANGULATION_H
#define _VTK_TRIANGULATION_H

// c++ includes
#include                  <algorithm>

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>

// VTK includes
#include                  <vtkCellArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkImageData.h>
#include                  <vtkPolyData.h>
#include                  <vtkUnstructuredGrid.h>

class vtkTriangulation : public Triangulation{

  public:
    
    vtkTriangulation();
    
    /// Specify the input VTK object representing a triangulation or a regular 
    /// grid.
    /// \param dataSet Input VTK object (vtkImageData, vtkPolyData or 
    /// vtkUnstructuredGrid).
    /// \return Returns 0 upon success, negative values otherwise.
    int setInputData(vtkDataSet *dataSet);
      
  protected:
  
    
  private:
    
};

#endif // _VTK_TRIANGULATION_H
