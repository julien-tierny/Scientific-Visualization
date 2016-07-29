/// \ingroup vtkWrappers
/// \class vtkReebSpace
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date October 2015.
///
/// \brief VTK wrapper for the reebSpace package.
///
/// VTK wrapping code for the @ReebSpace package.
/// \sa wtfit::ReebSpace
#ifndef _VTK_REEBSPACE_H
#define _VTK_REEBSPACE_H

// wtfit code includes
#include                  <ReebSpace.h>
#include                  <Wrapper.h>

// wtfit VTK helpers
#include                  <vtkOneSkeleton.h>
#include                  <vtkZeroSkeleton.h>
#include                  <vtkThreeSkeleton.h>

// VTK includes -- to adapt
#include                  <vtkCellArray.h>
#include                  <vtkCellData.h>
#include                  <vtkCharArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkInformationVector.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkUnstructuredGridAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a 
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK 
// class your wrapper should inherit.
class VTKFILTERSCORE_EXPORT vtkReebSpace 
  : public vtkUnstructuredGridAlgorithm, public Wrapper{

  public:
      
    static vtkReebSpace* New();
    
    vtkTypeMacro(vtkReebSpace, vtkUnstructuredGridAlgorithm);
    
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
    
    vtkSetMacro(Ucomponent, string);
    vtkGetMacro(Ucomponent, string);
    
    vtkSetMacro(Vcomponent, string);
    vtkGetMacro(Vcomponent, string);
    
    vtkGetMacro(VaryingConnectivity, bool);
    vtkSetMacro(VaryingConnectivity, bool);

    vtkGetMacro(VaryingValues, bool);
    vtkSetMacro(VaryingValues, bool);
    
    // 0-sheet options
    vtkGetMacro(ZeroSheetId, bool);
    vtkSetMacro(ZeroSheetId, bool);
    
    vtkGetMacro(ZeroSheetType, bool);
    vtkSetMacro(ZeroSheetType, bool);
    
    vtkGetMacro(ZeroSheetValue, bool);
    vtkSetMacro(ZeroSheetValue, bool);
    
    vtkGetMacro(ZeroSheetVertexId, bool);
    vtkSetMacro(ZeroSheetVertexId, bool);
    
    
    // 1-sheet options
    vtkGetMacro(OneSheetId, bool);
    vtkSetMacro(OneSheetId, bool);
    
    vtkGetMacro(OneSheetType, bool);
    vtkSetMacro(OneSheetType, bool);
    
    vtkGetMacro(OneSheetValue, bool);
    vtkSetMacro(OneSheetValue, bool);
    
    vtkGetMacro(OneSheetVertexId, bool);
    vtkSetMacro(OneSheetVertexId, bool);
    
    vtkGetMacro(OneSheetEdgeId, bool);
    vtkSetMacro(OneSheetEdgeId, bool);
    
    // 2-sheet options
    vtkGetMacro(TwoSheets, bool);
    vtkSetMacro(TwoSheets, bool);
    
    vtkGetMacro(TwoSheetCaseId, bool);
    vtkSetMacro(TwoSheetCaseId, bool);
    
    vtkGetMacro(TwoSheetValue, bool);
    vtkSetMacro(TwoSheetValue, bool);
    
    vtkGetMacro(TwoSheetParameterization, bool);
    vtkSetMacro(TwoSheetParameterization, bool);
    
    vtkGetMacro(TwoSheetId, bool);
    vtkSetMacro(TwoSheetId, bool);
    
    vtkGetMacro(TwoSheetEdgeId, bool);
    vtkSetMacro(TwoSheetEdgeId, bool);
    
    vtkGetMacro(TwoSheetTetId, bool);
    vtkSetMacro(TwoSheetTetId, bool);
    
    vtkGetMacro(TwoSheetEdgeType, bool);
    vtkSetMacro(TwoSheetEdgeType, bool);
    
    vtkGetMacro(ThreeSheetTetNumber, bool);
    vtkSetMacro(ThreeSheetTetNumber, bool);
    
    vtkGetMacro(ThreeSheetVertexNumber, bool);
    vtkSetMacro(ThreeSheetVertexNumber, bool);
    
    vtkGetMacro(ThreeSheetExpansion, bool);
    vtkSetMacro(ThreeSheetExpansion, bool);
    
    vtkGetMacro(ThreeSheetDomainVolume, bool);
    vtkSetMacro(ThreeSheetDomainVolume, bool);
    
    vtkGetMacro(ThreeSheetRangeArea, bool);
    vtkSetMacro(ThreeSheetRangeArea, bool);
    
    vtkGetMacro(ThreeSheetHyperVolume, bool);
    vtkSetMacro(ThreeSheetHyperVolume, bool);
    
    vtkGetMacro(SimplificationThreshold, double);
    vtkSetMacro(SimplificationThreshold, double);
    
    vtkGetMacro(SimplificationCriterion, int);
    vtkSetMacro(SimplificationCriterion, int);
    
  protected:
    
    vtkReebSpace();
    
    ~vtkReebSpace();
    
    virtual int FillOutputPortInformation(int port, vtkInformation *info){
      
      if(port == 0){
        // 0-sheets, corners of jacobi set segments
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      }
      else if(port == 1){
        // 1-sheets, jacobi sets
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      }
      else if(port == 2){
        // 2-sheets, fiber surfaces of jacobi sets
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      }
      else if(port == 3){
        // 3-sheets, pre-images of reeb space 2-sheets
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      }
      
      return 1;
    }
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **inputVector, vtkInformationVector *outputVector);
    
  private:
    
    bool                  UseAllCores;
    int                   ThreadNumber;
   
    bool                  VaryingConnectivity, VaryingValues;
    
    bool                  ZeroSheetValue, ZeroSheetVertexId, ZeroSheetType,
                          ZeroSheetId;
    bool                  OneSheetValue, OneSheetVertexId, OneSheetType,
                          OneSheetId, OneSheetEdgeId;
    bool                  TwoSheets,
                          TwoSheetCaseId, TwoSheetValue, 
                          TwoSheetParameterization,
                          TwoSheetId, TwoSheetEdgeId, TwoSheetTetId,
                          TwoSheetEdgeType;
    bool                  ThreeSheetVertexNumber, ThreeSheetTetNumber,
                          ThreeSheetExpansion,
                          ThreeSheetDomainVolume, ThreeSheetRangeArea,
                          ThreeSheetHyperVolume;

    bool                  perturbated_;
    int                   SimplificationCriterion;
    double                SimplificationThreshold;
    string                Ucomponent, Vcomponent;
    
    // connectivity preprocessing
    // jacobi set stuff
    vector<pair<int, int> > 
                          edgeList_;
    vector<vector<int> >  edgeStars_;
    vector<vector<pair<int, int> > >
                          edgeFanLinkEdgeList_;
    vector<vector<long long int> > 
                          edgeFans_;
    vector<vector<int> >  tetNeighbors_;
    vector<vector<int> >  vertexStars_;
    
    vector<int>           sosOffsets_;
    // reeb space traversal stuff
    vector<vector<int> >  vertexEdgeList_;
   
    // core data-structure
    ReebSpace             reebSpace_;
    
    // template base call 
    template <class dataTypeU, class dataTypeV> int baseCall(
      vtkUnstructuredGrid *input,
      vtkDataArray *uField,
      vtkDataArray *vField);
    
    template <class dataTypeU, class dataTypeV> int preProcess(
      vtkUnstructuredGrid *input,
      vtkDataArray *uField,
      vtkDataArray *vField);
    
    // base code features
    int doIt(vtkUnstructuredGrid *input, 
      vtkUnstructuredGrid *sheet0, vtkUnstructuredGrid *sheet1,
      vtkUnstructuredGrid *sheet2, vtkUnstructuredGrid *sheet3);
    
    bool needsToAbort();
    
    int updateProgress(const float &progress);
   
};

#endif // _VTK_REEBSPACE_H