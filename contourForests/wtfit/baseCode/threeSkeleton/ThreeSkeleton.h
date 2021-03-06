/// \ingroup baseCode
/// \class wtfit::ThreeSkeleton 
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date August 2015.
///
/// \brief %ThreeSkeleton processing package.
///
/// %ThreeSkeleton is a processing package that handles the 3-skeleton 
/// (tetrahedra) of a triangulation.
/// \sa Triangulation
/// \sa vtkTriangulation
/// \sa vtkThreeSkeleton

#ifndef _THREESKELETON_H
#define _THREESKELETON_H

// base code includes
#include                  <OneSkeleton.h>
#include                  <TwoSkeleton.h>
#include                  <ZeroSkeleton.h>
#include                  <Wrapper.h>

#include                  <algorithm>

namespace wtfit{
  
  class ThreeSkeleton : public Debug{

    public:
        
      ThreeSkeleton();
      
      ~ThreeSkeleton();
     
      /// Compute the list of edges of each cell of a triangulation.
      /// \param vertexNumber Number of vertices in the triangulation.
      /// \param cellNumber Number of maximum-dimensional cells in the 
      /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
      /// \param cellArray Pointer to a contiguous array of cells. Each entry 
      /// starts by the number of vertices in the cell, followed by the vertex
      /// identifiers of the cell.
      /// \param cellEdges Output edge lists. The size of this vector 
      /// will be equal to the number of cells in the mesh. Each entry will be
      /// a vector listing the edge identifiers of the entry's cell's 
      /// edges.
      /// \param edgeList Optional list of edges. If NULL, the function will 
      /// compute this list anyway and free the related memory upon return.
      /// If not NULL but pointing to an empty vector, the function will fill
      /// this empty vector (useful if this list needs to be used later on by 
      /// the calling program). If not NULL but pointing to a non-empty vector, 
      /// this function will use this vector as internal edge list. If this 
      /// vector is not empty but incorrect, the behavior is unspecified.
      /// \param vertexEdges Optional list of edges for each vertex. If NULL, 
      /// the function will compute this list anyway and free the related 
      /// memory upon return. If not NULL but pointing to an empty vector, 
      /// the function will fill this empty vector (useful if this list needs to
      /// be used later on by the calling program). If not NULL but pointing to 
      /// a non-empty vector, this function will use this vector as internal 
      /// vertex edge list. If this vector is not empty but incorrect, the 
      /// behavior is unspecified.
      /// \return Returns 0 upon success, negative values otherwise.
      int buildCellEdges(const int &vertexNumber,
        const int &cellNumber, 
        const long long int *cellArray,
        vector<vector<int> > &cellEdges,
        vector<pair<int, int> > *edgeList = NULL,
        vector<vector<int> > *vertexEdges = NULL) const ;
      
      /// Compute the list of cell-neighbors of each cell of a triangulation 
      /// (unspecified behavior if the input mesh is not a triangulation).
      /// This implementation is fast only if you already have the triangle 
      /// stars computed. Otherwise, please use 
      /// ThreeSkeleton::buildCellNeighborsFromVertices instead.
      /// \param vertexNumber Number of vertices in the triangulation.
      /// \param cellNumber Number of maximum-dimensional cells in the 
      /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
      /// \param cellArray Pointer to a contiguous array of cells. Each entry 
      /// starts by the number of vertices in the cell, followed by the vertex
      /// identifiers of the cell.
      /// \param cellNeighbors Output neighbor list. The size of this vector 
      /// will be equal to the number of cells in the mesh. Each entry will be
      /// a vector listing the cell identifiers of the entry's cell's 
      /// neighbors.
      /// \param triangleStars Optional list of triangle stars (list of 
      /// 3-dimensional cells connected to each triangle). If NULL, the 
      /// function will compute this list anyway and free the related memory
      /// upon return. If not NULL but pointing to an empty vector, the 
      /// function will fill this empty vector (useful if this list needs 
      /// to be used later on by the calling program). If not NULL but pointing
      /// to a non-empty vector, this function will use this vector as internal 
      /// triangle star list. If this vector is not empty but incorrect, the 
      /// behavior is unspecified.
      /// \return Returns 0 upon success, negative values otherwise.
      int buildCellNeighborsFromTriangles(const int &vertexNumber, 
        const int &cellNumber,
        const long long int *cellArray,
        vector<vector<int> > &cellNeighbors,
        vector<vector<int> > *triangleStars = NULL) const;

      /// Compute the list of cell-neighbors of each cell of a triangulation 
      /// (unspecified behavior if the input mesh is not a triangulation).
      /// \param vertexNumber Number of vertices in the triangulation.
      /// \param cellNumber Number of maximum-dimensional cells in the 
      /// triangulation (number of tetrahedra in 3D, triangles in 2D, etc.)
      /// \param cellArray Pointer to a contiguous array of cells. Each entry 
      /// starts by the number of vertices in the cell, followed by the vertex
      /// identifiers of the cell.
      /// \param cellNeighbors Output neighbor list. The size of this vector 
      /// will be equal to the number of cells in the mesh. Each entry will be
      /// a vector listing the cell identifiers of the entry's cell's 
      /// neighbors.
      /// \param vertexStars Optional list of vertex stars (list of 
      /// 3-dimensional cells connected to each vertex). If NULL, the 
      /// function will compute this list anyway and free the related memory
      /// upon return. If not NULL but pointing to an empty vector, the 
      /// function will fill this empty vector (useful if this list needs 
      /// to be used later on by the calling program). If not NULL but pointing
      /// to a non-empty vector, this function will use this vector as internal 
      /// vertex star list. If this vector is not empty but incorrect, the 
      /// behavior is unspecified.
      /// \return Returns 0 upon success, negative values otherwise.
      int buildCellNeighborsFromVertices(const int &vertexNumber, 
        const int &cellNumber,
        const long long int *cellArray,
        vector<vector<int> > &cellNeighbors,
        vector<vector<int> > *vertexStars = NULL) const;
        
    protected:
    
  };
}

// if the package is not a template, comment the following line
// #include                  <ThreeSkeleton.cpp>

#endif // THREESKELETON_H
