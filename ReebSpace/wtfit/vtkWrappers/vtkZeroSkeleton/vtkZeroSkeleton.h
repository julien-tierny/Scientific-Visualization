/// \ingroup vtkWrappers
/// \class vtkZeroSkeleton
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief VTK helper that processes the 0-skeleton (vertices) of a data-set.
/// \sa vtkTriangulation
/// \sa wtfit::Triangulation
/// \sa wtfit::ZeroSkeleton

#ifndef _VTK_ZERO_SKELETON_H
#define _VTK_ZERO_SKELETON_H

// c++ includes
#include                  <algorithm>

// base code includes
#include                  <Debug.h>
#include                  <ZeroSkeleton.h>
#include                  <vtkOneSkeleton.h>

// VTK includes
#include                  <vtkCellArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkGenericCell.h>
#include                  <vtkIdList.h>
#include                  <vtkPolyData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkUnstructuredGrid.h>

// TODO:
// make a VTK filter if needed

class vtkZeroSkeleton : public Wrapper{

  public:
    
    vtkZeroSkeleton();
   
    /// Compute the list of edges connected to each vertex of a vtkDataSet.
    /// \param input Input data-set.
    /// \param vertexEdges Output vertex links. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the identifiers of the edges connected to the entry's
    /// vertex.
    /// \param edgeList Optional list of edges. If NULL, the function will 
    /// compute this list anyway and free the related memory upon return.
    /// If not NULL but pointing to an empty vector, the function will fill this
    /// empty vector (useful if this list needs to be used later on by the 
    /// calling program). If not NULL but pointing to a non-empty vector, this
    /// function will use this vector as internal edge list. If this vector is
    /// not empty but incorrect, the behavior is unspecified.
    /// \param isTriangulation Optional flag that speeds up computation if the 
    /// input mesh is indeed a valid triangulation (unspecified behavior 
    /// otherwise).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexEdges(vtkDataSet *input,
      vector<vector<int> > &vertexEdges,
      vector<pair<int, int> > *edgeList = NULL,
      const bool &isTriangulation = false) const;
    
    /// Compute the link of each vertex of a vtkDataSet.
    /// \param input Input data-set.
    /// \param vertexLinks Output vertex links. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the simplices of the link of the entry's vertex. In 
    /// particular, this vector contains, for each simplex, the number of 
    /// vertices in the simplex (triangles: 3, edges: 2) followed by the 
    /// corresponding vertex identifiers.
    /// \param vertexStars Optional list of vertex stars (list of 
    /// 3-dimensional cells connected to each vertex). If NULL, the 
    /// function will compute this list anyway and free the related memory
    /// upon return. If not NULL but pointing to an empty vector, the 
    /// function will fill this empty vector (useful if this list needs 
    /// to be used later on by the calling program). If not NULL but pointing
    /// to a non-empty vector, this function will use this vector as internal 
    /// vertex star list. If this vector is not empty but incorrect, the 
    /// behavior is unspecified.
    /// \param isTriangulation Optional flag that speeds up computation if the 
    /// input mesh is indeed a valid triangulation (unspecified behavior 
    /// otherwise).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexLinks(vtkDataSet *input, 
      vector<vector<long long int > > &vertexLinks,
      vector<vector<int> > *vertexStars = NULL,
      const bool &isTriangulation = false) const;
    
    /// Compute the list of neighbors of each vertex of a vtkDataSet.
    /// \param input Input data-set.
    /// \param vertexNeighbors Output neighbor list. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the vertex identifiers of the entry's vertex' 
    /// neighbors.
    /// \param isTriangulation Optional flag that speeds up computation if the 
    /// input mesh is indeed a valid triangulation (unspecified behavior 
    /// otherwise).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexNeighbors(vtkDataSet *input, 
      vector<vector<int> > &vertexNeighbors,
      const bool &isTriangulation = false) const;

    /// Compute the star of each vertex of a vtkDataSet.
    /// \param input Input data-set.
    /// \param vertexStars Output vertex stars. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the identifiers of the maximum-dimensional cells 
    /// (3D: tetrahedra, 2D: triangles, etc.) connected to the entry's vertex.
    /// \param isTriangulation Optional flag that speeds up computation if the 
    /// input mesh is indeed a valid triangulation (unspecified behavior 
    /// otherwise).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildVertexStars(vtkDataSet *input, 
      vector<vector<int> > &vertexStars, 
      const bool &isTriangulation = false) const;
     
    /// Compute the link of each vertex of a triangulation represented by a
    /// vtkPolyData object (unspecified behavior if the input mesh is
    /// not a valid triangulation).
    /// \param input Input data-set.
    /// \param vertexLinks Output vertex links. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the simplices of the link of the entry's vertex. In 
    /// particular, this vector contains, for each simplex, the number of 
    /// vertices in the simplex (triangles: 3, edges: 2) followed by the 
    /// corresponding vertex identifiers.
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
    int buildTriangulationVertexLinks(vtkPolyData *input, 
      vector<vector<long long int > > &vertexLinks,
      vector<vector<int> > *vertexStars = NULL) const{
      
      ZeroSkeleton zeroSkeleton;
      zeroSkeleton.setWrapper(this);
      return zeroSkeleton.buildVertexLinks(input->GetNumberOfPoints(),
        input->GetNumberOfCells(), input->GetPolys()->GetPointer(),
        vertexLinks, vertexStars);
    }
    
    /// Compute the link of each vertex of a triangulation represented by a
    /// vtkUnstructuredGrid object (unspecified behavior if the input mesh is
    /// not a valid triangulation).
    /// \param input Input data-set.
    /// \param vertexLinks Output vertex links. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the simplices of the link of the entry's vertex. In 
    /// particular, this vector contains, for each simplex, the number of 
    /// vertices in the simplex (triangles: 3, edges: 2) followed by the 
    /// corresponding vertex identifiers.
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
    int buildTriangulationVertexLinks(vtkUnstructuredGrid *input, 
      vector<vector<long long int > > &vertexLinks,
      vector<vector<int> > *vertexStars = NULL) const{
      
      ZeroSkeleton zeroSkeleton;
      zeroSkeleton.setWrapper(this);
      return zeroSkeleton.buildVertexLinks(input->GetNumberOfPoints(),
        input->GetNumberOfCells(), input->GetCells()->GetPointer(),
        vertexLinks, vertexStars);
    }
      
    /// Compute the list of neighbors of each vertex of a vtkPolyData. 
    /// Unspecified behavior if the input mesh is not a valid triangulation).
    /// \param input Input data-set.
    /// \param vertexNeighbors Output neighbor list. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the vertex identifiers of the entry's vertex' 
    /// neighbors.
    /// \param edgeList Optional list of edges. If NULL, the function will 
    /// compute this list anyway and free the related memory upon return.
    /// If not NULL but pointing to an empty vector, the function will fill this
    /// empty vector (useful if this list needs to be used later on by the 
    /// calling program). If not NULL but pointing to a non-empty vector, this
    /// function will use this vector as internal edge list. If this vector is
    /// not empty but incorrect, the behavior is unspecified      
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangulationVertexNeighbors(vtkPolyData *input,
      vector<vector<int> > &vertexNeighbors,
      vector<pair<int, int> > *edgeList = NULL) const{
    
      ZeroSkeleton skeleton;
      skeleton.setWrapper(this);
      return skeleton.buildVertexNeighbors(input->GetNumberOfPoints(), 
        input->GetNumberOfCells(), input->GetPolys()->GetPointer(),
        vertexNeighbors, edgeList);
    }
    
    /// Compute the list of neighbors of each vertex of a vtkUnstructuredGrid. 
    /// Unspecified behavior if the input mesh is not a valid triangulation).
    /// \param input Input data-set.
    /// \param vertexNeighbors Output neighbor list. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be 
    /// a vector listing the vertex identifiers of the entry's vertex' 
    /// neighbors.
    /// \param edgeList Optional list of edges. If NULL, the function will 
    /// compute this list anyway and free the related memory upon return.
    /// If not NULL but pointing to an empty vector, the function will fill this
    /// empty vector (useful if this list needs to be used later on by the 
    /// calling program). If not NULL but pointing to a non-empty vector, this
    /// function will use this vector as internal edge list. If this vector is
    /// not empty but incorrect, the behavior is unspecified      
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangulationVertexNeighbors(vtkUnstructuredGrid *input,
      vector<vector<int> > &vertexNeighbors,
      vector<pair<int, int> > *edgeList = NULL) const{
    
      ZeroSkeleton skeleton;
      skeleton.setWrapper(this);
      return skeleton.buildVertexNeighbors(input->GetNumberOfPoints(), 
        input->GetNumberOfCells(), input->GetCells()->GetPointer(),
        vertexNeighbors, edgeList);
    }
    
    /// Compute the star of each vertex of a triangulation represented by a 
    /// vtkPolyData object. Unspecified behavior if the input mesh is not a 
    /// valid triangulation.
    /// \param input Input data-set.
    /// \param vertexStars Output vertex stars. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the identifiers of the maximum-dimensional cells 
    /// (3D: tetrahedra, 2D: triangles, etc.) connected to the entry's vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangulationVertexStars(vtkPolyData *input,
      vector<vector<int> > &vertexStars) const{
        
      ZeroSkeleton zeroSkeleton;
      return zeroSkeleton.buildVertexStars(input->GetNumberOfPoints(),
        input->GetNumberOfCells(),
        input->GetPolys()->GetPointer(),
        vertexStars);
    }
    
    /// Compute the star of each vertex of a triangulation represented by a 
    /// vtkUnstructuredGrid object. Unspecified behavior if the input mesh is 
    /// not a valid triangulation.
    /// \param input Input data-set.
    /// \param vertexStars Output vertex stars. The size of this vector 
    /// will be equal to the number of vertices in the mesh. Each entry will be
    /// a vector listing the identifiers of the maximum-dimensional cells 
    /// (3D: tetrahedra, 2D: triangles, etc.) connected to the entry's vertex.
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangulationVertexStars(vtkUnstructuredGrid *input,
      vector<vector<int> > &vertexStars) const{
        
      ZeroSkeleton zeroSkeleton;
      return zeroSkeleton.buildVertexStars(input->GetNumberOfPoints(),
        input->GetNumberOfCells(),
        input->GetCells()->GetPointer(),
        vertexStars);
    }
      
    
      
  protected:
  
    
  private:
    
    // empty wrapping to VTK for now
    bool needsToAbort(){ return false;};
    
    int updateProgress(const float &progress) {return 0;};
};

#endif // _VTK_ZERO_SKELETON_H
