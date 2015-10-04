/* 
 * file:                  vtkFiberSurfaceFilter.h
 * description:           VTK wrapper for the FiberSurface package.
 * author(s):             Pavol Klacansky <pavol@klacansky.com>.
 *                        Julien Tierny <julien.tierny@lip6.fr>
 * date:                  March 2015.
 */

/// \ingroup vtkWrappers
/// \class vtkFiberSurfaceFilter
/// \brief VTK wrapper for the FiberSurface package.
///
/// VTK wrapping code for the @FiberSurface package.
/// \sa wtfit::FiberSurface
#ifndef _VTK_FIBER_SURFACE_H
#define _VTK_FIBER_SURFACE_H

#include                  <cfloat>
#include                  <cinttypes>
#include                  <cmath>
#include                  <cstdint>
#include                  <iostream>
#include                  <string>

#include                  <RangeDrivenOctree.h>

// ParaView include
#include                  <vtkCleanUnstructuredGrid.h>

// VTK includes 
#include                  <vtkAppendPolyData.h>
#include                  <vtkCell.h>
#include                  <vtkCellArray.h>
#include                  <vtkCellData.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataObject.h>
#include                  <vtkDataSetSurfaceFilter.h>
#include                  <vtkDemandDrivenPipeline.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkPointData.h>
#include                  <vtkPoints.h>
#include                  <vtkPolyDataAlgorithm.h>
#include                  <vtkSmartPointer.h>
#include                  <vtkUnstructuredGrid.h>
#include                  <vtkXMLUnstructuredGridWriter.h>

extern "C" {
    #include "bvh.h"
}

using namespace std;

class VTKFILTERSCORE_EXPORT vtkFiberSurfaceFilter 
  : public vtkPolyDataAlgorithm, public Wrapper{

  public:
      
    static vtkFiberSurfaceFilter* New();
    
    vtkTypeMacro(vtkFiberSurfaceFilter, vtkPolyDataAlgorithm);
   
    vtkSetMacro(DataUComponent, string);
    vtkGetMacro(DataUComponent, string);

    vtkSetMacro(DataVComponent, string);
    vtkGetMacro(DataVComponent, string);
    
    vtkSetMacro(PolygonUComponent, string);
    vtkGetMacro(PolygonUComponent, string);

    vtkSetMacro(PolygonVComponent, string);
    vtkGetMacro(PolygonVComponent, string); 

    vtkSetMacro(ThreadNumber, int);
    vtkGetMacro(ThreadNumber, int);    
    
    enum vtkFiberSurfaceImplementation {
      REGULAR, OCTREE, BVH};
    vtkSetMacro(Implementation, int);
    
    enum vtkFiberSurfaceThreadingStrategy{
      AUTOMATIC, DOMAIN_TETRAHEDRA, POLYGON_EDGES
    };
    vtkSetMacro(ThreadStrategy, int);
    vtkSetMacro(ThreadBalance, float);
   
    enum vtkFiberSurfaceAlgorithm {
      SIMPLE, GREY, GREY_MINIMAL
    };
    vtkSetMacro(Algorithm, int);  
    
    vtkSetMacro(Manifold, bool);
    vtkSetMacro(PolygonSegmentation, bool);
    vtkSetMacro(ShowCandidates, bool);
    vtkSetMacro(VisibleFibers, bool);
    
    vtkSetMacro(debugLevel_, int);
    
    void SetBvhMinimumTetNumber(int number){
      BvhMinimumTetNumber = number;
      if(BVH_){
        bvh_free(BVH_);
        BVH_ = NULL;
      }
      Modified();
    }
    
    void SetOctreeRangeRatio(float ratio){
      octreeMinimumRangeAreaRatio_ = ratio;
      if(octree_){
        switch(dataType0_){
          case VTK_CHAR:
            switch(dataType1_){
              vtkTemplateMacro((
                delete ((RangeDrivenOctree<char, VTK_TT> *) octree_)
              ));
            }
            break;
          case VTK_DOUBLE:
            switch(dataType1_){
              vtkTemplateMacro((
                delete ((RangeDrivenOctree<double, VTK_TT> *) octree_)
              ));
            }
            break;
          case VTK_FLOAT:
            switch(dataType1_){
              vtkTemplateMacro((
                delete ((RangeDrivenOctree<float, VTK_TT> *) octree_)
              ));
            }
            break;
          case VTK_INT:
            switch(dataType1_){
              vtkTemplateMacro((
                delete ((RangeDrivenOctree<int, VTK_TT> *) octree_)
              ));
            }
            break;
          case VTK_UNSIGNED_CHAR:
            switch(dataType1_){
              vtkTemplateMacro((
                delete ((RangeDrivenOctree<unsigned char, VTK_TT> *) octree_)
              ));
            }
            break;
        }
        octree_ = NULL;
      }
      Modified();
    }    

    bool needsToAbort(){
      return GetAbortExecute();
    }
    
    int updateProgress(const float &progress){
      UpdateProgress(progress);
      return 0;
    }
    
  protected:
    
    template<class dataType0, class dataType1> int buildBVH(
      vtkUnstructuredGrid *input, void *u, void *v);
    
    template<class dataType0, class dataType1> int buildOctree(
      vtkUnstructuredGrid *input, void *u, void *v);
    
    double doIt(vtkUnstructuredGrid *input, vtkUnstructuredGrid *polygon, 
      vtkPolyData *output);

    double doItBVH(vtkUnstructuredGrid *input, 
      vtkUnstructuredGrid *polygon, vtkPolyData *output);
    
    double doItOctree(vtkUnstructuredGrid *input, 
      vtkUnstructuredGrid *polygon, 
      vtkPolyData *output);
   
    int reduce();
    
    vtkFiberSurfaceFilter();
    
    ~vtkFiberSurfaceFilter();
    
    int FillInputPortInformation(int port, vtkInformation *info);
    
    int RequestData(vtkInformation *request, 
      vtkInformationVector **in_vec, vtkInformationVector *out_vec);
    
    
  private:
    
    bool                  Manifold, PolygonSegmentation,
                          ShowCandidates, VisibleFibers;
    float                 ThreadBalance, octreeMinimumRangeAreaRatio_;
    int                   BvhMinimumTetNumber;
    int                   ThreadNumber, ThreadStrategy;
    int                   Implementation;
    int                   Algorithm;
    int                   dataType0_, dataType1_;
    void                  *octree_;
    vtkUnstructuredGrid   *input;
    vtkPolyData           *output;
    vector<vtkSmartPointer<vtkPolyData> > 
                          threadedOutput_;
    vector<vtkSmartPointer<vtkPoints> >   
                          threadedPoints_;
    vector<vtkSmartPointer<vtkCellArray> > 
                          threadedCells_;
    vector<vtkSmartPointer<vtkIntArray> >
                          threadedEdgeIds_;
    vector<vtkSmartPointer<vtkFloatArray> >
                          threadedTextureCoordinates_;
    vtkSmartPointer<vtkAppendPolyData>     
                          appender_;
    struct Bvh            *BVH_;
    
    static const char name[];
   
    string                DataUComponent, DataVComponent;
    string                PolygonUComponent, PolygonVComponent;
    
};

#endif // _VTK_FIBER_SURFACE_H
