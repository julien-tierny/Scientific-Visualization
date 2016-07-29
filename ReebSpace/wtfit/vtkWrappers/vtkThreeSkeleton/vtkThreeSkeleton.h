/// \ingroup vtkWrappers
/// \class vtkThreeSkeleton
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date August 2015.
///
/// \brief VTK helper that processes the 3-skeleton (tetrahedra) of a data-set.
///
/// VTK wrapping code for the @ThreeSkeleton package.
/// \sa vtkTriangulation
/// \sa wtfit::Triangulation
/// \sa wtfit::ThreeSkeleton
#ifndef _VTK_THREESKELETON_H
#define _VTK_THREESKELETON_H

// wtfit code includes
#include                  <ThreeSkeleton.h>
#include                  <Wrapper.h>
#include                  <vtkZeroSkeleton.h>
#include                  <vtkOneSkeleton.h>

// VTK includes -- to adapt
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkGenericCell.h>
#include                  <vtkInformation.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkSmartPointer.h>


class VTKFILTERSCORE_EXPORT vtkThreeSkeleton 
  : public vtkDataSetAlgorithm, public Wrapper{

  public:
      
    static vtkThreeSkeleton* New();
    
    vtkTypeMacro(vtkThreeSkeleton, vtkDataSetAlgorithm);
    
    // default wtfit setters
    vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = OsCall::getNumberOfCores();
      }
      Modified();
    }
    
    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }   
    
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default wtfit setters
   
    /// Compute the list of edges of each cell of a vtkDataSet.
    /// \param input Input data-set.
    /// \param cellEdges Output edge lists. The size of this vector 
    /// will be equal to the number of cells in the mesh. Each entry will be
    /// a vector listing the edge identifiers of the entry's cell's 
    /// edges.
    /// \param edgeList Optional list of edges. If NULL, the function will 
    /// compute this list anyway and free the related memory upon return.
    /// If not NULL but pointing to an empty vector, the function will fill this
    /// empty vector (useful if this list needs to be used later on by the 
    /// calling program). If not NULL but pointing to a non-empty vector, this
    /// function will use this vector as internal edge list. If this vector is
    /// not empty but incorrect, the behavior is unspecified.
    /// \param vertexEdges Optional list of edges for each vertex. If NULL, the
    /// function will compute this list anyway and free the related memory upon
    /// return. If not NULL but pointing to an empty vector, the function will
    /// fill this empty vector (useful if this list needs to be used later on by
    /// the calling program). If not NULL but pointing to a non-empty vector, 
    /// this function will use this vector as internal vertex edge list. If 
    /// this vector is not empty but incorrect, the behavior is unspecified.
    /// \param isTriangulation Optional flag that speeds up computation if the 
    /// input mesh is indeed a valid triangulation (unspecified behavior 
    /// otherwise).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildCellEdges(vtkDataSet *input,
      vector<vector<int> > &cellEdges,
      vector<pair<int, int> > *edgeList = NULL,
      vector<vector<int> > *vertexEdges = NULL,
      const bool &isTriangulation = false) const;
        
    /// Compute the list of neighbors of each cell of a vtkDataSet.
    /// \param input Input data-set.
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
    /// \param triangleStars Optional list of triangle stars (list of 
    /// 3-dimensional cells connected to each triangle). If NULL, the 
    /// function will compute this list anyway and free the related memory
    /// upon return. If not NULL but pointing to an empty vector, the 
    /// function will fill this empty vector (useful if this list needs 
    /// to be used later on by the calling program). If not NULL but pointing
    /// to a non-empty vector, this function will use this vector as internal 
    /// triangle star list. If this vector is not empty but incorrect, the 
    /// behavior is unspecified.
    /// \param isTriangulation Optional flag that speeds up computation if the 
    /// input mesh is indeed a valid triangulation (unspecified behavior 
    /// otherwise).
    /// \return Returns 0 upon success, negative values otherwise.
    int buildCellNeighbors(vtkDataSet *input,
      vector<vector<int> > &cellNeighbors, 
      vector<vector<int> > *vertexStars = NULL,
      vector<vector<int> > *triangleStars = NULL,
      const bool &isTriangulation = false) const;

    /// Compute the list of neighbors of each cell of a triangulation 
    /// represented by.a vtkPolyData object (unspecified behavior if the input 
    /// mesh is not a triangulation).
    /// \param input Input vtkPolyData object.
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
    int buildTriangulationCellNeighbors(vtkPolyData *input,
      vector<vector<int> > &cellNeighbors, 
      vector<vector<int> > *vertexStars = NULL,
      vector<vector<int> > *triangleStars = NULL) const{
       
      ThreeSkeleton threeSkeleton;
      
      threeSkeleton.setWrapper(this);
      
      if((triangleStars)&&(!vertexStars)){
        return threeSkeleton.buildCellNeighborsFromTriangles(
          input->GetNumberOfPoints(),
          input->GetNumberOfCells(),
          input->GetPolys()->GetPointer(),
          cellNeighbors, triangleStars);
      }
      else{ 
        return threeSkeleton.buildCellNeighborsFromVertices(
          input->GetNumberOfPoints(),
          input->GetNumberOfCells(),
          input->GetPolys()->GetPointer(),
          cellNeighbors, vertexStars);
      }
    }
    
    /// Compute the list of neighbors of each cell of a triangulation 
    /// represented by.a vtkUnstructuredGrid object (unspecified behavior if
    /// the input mesh is not a triangulation).
    /// \param input Input vtkUnstructuredGrid object.
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
    int buildTriangulationCellNeighbors(vtkUnstructuredGrid *input,
      vector<vector<int> > &cellNeighbors, 
      vector<vector<int> > *vertexStars = NULL,
      vector<vector<int> > *triangleStars = NULL) const{
       
      ThreeSkeleton threeSkeleton;
      
      threeSkeleton.setWrapper(this);
      
      if((triangleStars)&&(!vertexStars)){
        return threeSkeleton.buildCellNeighborsFromTriangles(
          input->GetNumberOfPoints(),
          input->GetNumberOfCells(),
          input->GetCells()->GetPointer(),
          cellNeighbors, triangleStars);
      }
      else{ 
        return threeSkeleton.buildCellNeighborsFromVertices(
          input->GetNumberOfPoints(),
          input->GetNumberOfCells(),
          input->GetCells()->GetPointer(),
          cellNeighbors, vertexStars);
      }
    }
    
  protected:
    
    vtkThreeSkeleton();
    
    ~vtkThreeSkeleton();
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);
    
    
  private:
    
    bool                  UseAllCores;
    int                   ThreadNumber;
    
    // base code features
    virtual int doIt(vtkDataSet *input, vtkDataSet *output);
    
    bool needsToAbort();
    
    int updateProgress(const float &progress);
   
};

#endif // _VTK_THREESKELETON_H