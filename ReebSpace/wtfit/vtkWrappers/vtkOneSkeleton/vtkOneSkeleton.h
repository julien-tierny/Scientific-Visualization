/// \ingroup vtkWrappers
/// \class vtkOneSkeleton
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date November 2014.
///
/// \brief VTK helper that processes the 1-skeleton (edges) of a data-set.
/// \sa vtkTriangulation
/// \sa wtfit::Triangulation
/// \sa wtfit::OneSkeleton

#ifndef _VTK_ONE_SKELETON_H
#define _VTK_ONE_SKELETON_H

// c++ includes
#include                  <algorithm>

// base code includes
#include                  <Debug.h>
#include                  <OneSkeleton.h>

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

class vtkOneSkeleton : public Wrapper{

  public:
    
    vtkOneSkeleton();
   
    /// Compute the list of edges of a vtkDataSet.
    /// \param input Input data-set.
    /// \param edgeList Output edge list (each entry is an ordered pair of 
    /// vertex identifiers).
    /// \param isTriangulation Optional flag that speeds up computation if the 
    /// input mesh is indeed a valid triangulation (unspecified behavior 
    /// otherwise).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeList(vtkDataSet *input, 
      vector<pair<int, int> > &edgeList, 
      const bool &isTriangulation = false) const;
      
    /// Compute the lists of edges of multiple triangulations (unspecified 
    /// behavior if the input meshes are not valid triangulations).
    /// \param input Input vtkUnstructuredGrid object.
    /// \param edgeList Output edge list (each entry is an ordered pair of 
    /// vertex identifiers).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildEdgeLists(const vector<vector<long long int> > &cellArays,
      vector<vector<pair<int, int> > > &edgeLists) const{
      
      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeLists(cellArays, edgeLists);
    }
    
    /// Compute the 3-star of all the edges of a vtkDataSet (for each edge, 
    /// list of the 3-dimensional cells connected to it).
    /// \param input Input data-set.
    /// \param starList Output list of 3-stars. The size of this vector will
    /// be equal to the number of edges in the mesh. Each entry stores a vector
    /// that lists the identifiers of all 3-dimensional cells connected to the 
    /// entry's edge.
    /// \param edgeList Optional list of edges. If NULL, the function will 
    /// compute this list anyway and free the related memory upon return.
    /// If not NULL but pointing to an empty vector, the function will fill this
    /// empty vector (useful if this list needs to be used later on by the 
    /// calling program). If not NULL but pointing to a non-empty vector, this
    /// function will use this vector as internal edge list. If this vector is
    /// not empty but incorrect, the behavior is unspecified.
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
    int buildEdgeStars(vtkDataSet *input,
      vector<vector<int> > &starList,
      vector<pair<int, int> > *edgeList = NULL,
      vector<vector<int> > *vertexStars = NULL,
      const bool &isTriangulation = false) const;
     
    /// Compute the list of edges of a triangulation represented by a 
    /// vtkPolyData object (unspecified behavior if the input mesh is not a 
    /// triangulation).
    /// \param input Input vtkPolyData object.
    /// \param edgeList Output edge list (each entry is an ordered pair of 
    /// vertex identifiers).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangulationEdgeList(vtkPolyData *input,
      vector<pair<int, int> > &edgeList) const{
      
      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeList(input->GetNumberOfPoints(), 
        input->GetNumberOfCells(), input->GetPolys()->GetPointer(), 
        edgeList);
        
      return 0;
    }
    
    /// Compute the list of edges of a triangulation represented by a 
    /// vtkUnstructuredGrid object (unspecified behavior if the input mesh is
    /// not a triangulation).
    /// \param input Input vtkUnstructuredGrid object.
    /// \param edgeList Output edge list (each entry is an ordered pair of 
    /// vertex identifiers).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildTriangulationEdgeList(vtkUnstructuredGrid *input,
      vector<pair<int, int> > &edgeList) const{
      
      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeList(input->GetNumberOfPoints(), 
        input->GetNumberOfCells(), input->GetCells()->GetPointer(), 
        edgeList);
        
      return 0;
    }
    
    /// Compute the 3-star of all the edges of a vtkPolyData (for each edge, 
    /// list of the 3-dimensional cells connected to it). Unspecified behavior 
    /// if the input mesh is not a valid triangulation).
    /// \param input Input data-set.
    /// \param starList Output list of 3-stars. The size of this vector will
    /// be equal to the number of edges in the mesh. Each entry stores a vector
    /// that lists the identifiers of all 3-dimensional cells connected to the 
    /// entry's edge.
    /// \param edgeList Optional list of edges. If NULL, the function will 
    /// compute this list anyway and free the related memory upon return.
    /// If not NULL but pointing to an empty vector, the function will fill this
    /// empty vector (useful if this list needs to be used later on by the 
    /// calling program). If not NULL but pointing to a non-empty vector, this
    /// function will use this vector as internal edge list. If this vector is
    /// not empty but incorrect, the behavior is unspecified.
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
    int buildTriangulationEdgeStars(vtkPolyData *input,
      vector<vector<int> > &starList,
      vector<pair<int, int> > *edgeList = NULL,
      vector<vector<int> > *vertexStars = NULL) const{
      
      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeStars(input->GetNumberOfPoints(),
        input->GetNumberOfCells(), input->GetPolys()->GetPointer(),
        starList, edgeList, vertexStars);
    }
    
    /// Compute the 3-star of all the edges of a vtkUnstructuredGrid (for each
    /// edge, list of the 3-dimensional cells connected to it). Unspecified 
    /// behavior if the input mesh is not a valid triangulation).
    /// \param input Input data-set.
    /// \param starList Output list of 3-stars. The size of this vector will
    /// be equal to the number of edges in the mesh. Each entry stores a vector
    /// that lists the identifiers of all 3-dimensional cells connected to the 
    /// entry's edge.
    /// \param edgeList Optional list of edges. If NULL, the function will 
    /// compute this list anyway and free the related memory upon return.
    /// If not NULL but pointing to an empty vector, the function will fill this
    /// empty vector (useful if this list needs to be used later on by the 
    /// calling program). If not NULL but pointing to a non-empty vector, this
    /// function will use this vector as internal edge list. If this vector is
    /// not empty but incorrect, the behavior is unspecified.
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
    int buildTriangulationEdgeStars(vtkUnstructuredGrid *input,
      vector<vector<int> > &starList,
      vector<pair<int, int> > *edgeList = NULL,
      vector<vector<int> > *vertexStars = NULL) const{
      
      OneSkeleton oneSkeleton;
      oneSkeleton.setWrapper(this);
      return oneSkeleton.buildEdgeStars(input->GetNumberOfPoints(),
        input->GetNumberOfCells(), input->GetCells()->GetPointer(),
        starList, edgeList, vertexStars);
    }
    
      
  protected:
  
    
  private:
    
    // empty wrapping to VTK for now
    bool needsToAbort(){ return false;};
    
    int updateProgress(const float &progress) {return 0;};
};

#endif // _VTK_ONE_SKELETON_H
