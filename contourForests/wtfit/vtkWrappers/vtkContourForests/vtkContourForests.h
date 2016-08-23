/// \ingroup vtkWrappers
/// \class vtkContourForests
/// \author Guillaume Favelier <guillaume.favelier@gmail.com>.
/// \date February 2016
///
/// \brief VTK wrapper for the contourTree package.
///
/// VTK wrapper for the @ContourTree package.
/// \sa wtfit::ContourTree

#ifndef _VTK_CONTOURTREE_H
#define _VTK_CONTOURTREE_H

// base code includes
#include <ParallelContourTree.h>
#include <ContourTree.h>
#include<Geometry.h>
#include <Wrapper.h>

// vtk wrapper includes
#include <vtkTriangulation.h>

// VTK includes
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkConnectivityFilter.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTable.h>

struct _persistenceCmp2 {
   _persistenceCmp2(const vector<double>* vertexScalars)
   {
      vertexScalars_ = vertexScalars;
   }

   bool operator()(const pair<pair<int, int>, double>& p0, const pair<pair<int, int>, double>& p1)
   {
      return (*vertexScalars_)[p0.first.first] < (*vertexScalars_)[p1.first.first];
   }

   const vector<double>* vertexScalars_;
};

enum class TreeComponent { ARC = -1, LOCAL_MINIMUM, SADDLE1, SADDLE2, LOCAL_MAXIMUM };

enum class TreeType { MERGE_TREE = 0, SPLIT_TREE = 1, CONTOUR_TREE = 2 };

enum class ArcType { MIN_ARC = 0, MAX_ARC, SADDLE1_ARC, SADDLE2_ARC, SADDLE1_SADDLE2_ARC };

enum class NodeType { LOCAL_MINIMUM = 0, SADDLE1, SADDLE2, LOCAL_MAXIMUM, REGULAR };

class VTKFILTERSCORE_EXPORT vtkContourForests : public vtkDataSetAlgorithm, public Wrapper
{
  public:
   static vtkContourForests* New();

   vtkTypeMacro(vtkContourForests, vtkDataSetAlgorithm);

   // default wtfit setters
   vtkSetMacro(debugLevel_, int);

   vtkSetMacro(FieldId, int);

   void SetThreads();
   void SetThreadNumber(int threadNumber);
   void SetDebugLevel(int d);
   void SetUseAllCores(bool onOff);
   // end of default wtfit setters

   vtkGetMacro(scalarField_, string);
   void SetScalarField(string scalarField);

   vtkGetMacro(useInputOffsetScalarField_, int);
   void SetUseInputOffsetScalarField(bool onOff);

   vtkSetMacro(inputOffsetScalarFieldName_, string);
   vtkGetMacro(inputOffsetScalarFieldName_, string);

   void CheckConformance(int checkConformance);
   void SetTreeType(int tree);

   void EnableVaryingMeshGeometry(bool state);
   void EnableVaryingMeshConnectivity(bool state);
   void EnableVaryingDataValues(bool state);

   void ShowMin(bool state);
   void ShowMax(bool state);
   void ShowSaddle1(bool state);
   void ShowSaddle2(bool state);
   void ShowSegmentation(bool segm);

   void ShowArc(bool state);
   void SetArcResolution(int arcResolution);
   void SetPartitionNumber(int partitionNum);
   void SetLessPartition(bool l);

   void SetSkeletonSmoothing(double skeletonSmooth);

   void ShowPersistenceCurve(bool state);
   void ShowPersistenceDiagram(bool state);

   void SetSimplificationType(int type);

   void SetSimplificationThreshold(double simplificationThreshold);

   void SetSimplifyPersistenceDiagram(bool state);

  protected:
   vtkContourForests();
   ~vtkContourForests();

   /// VTK Interface ///
   int RequestData(vtkInformation* request, vtkInformationVector** inputVector,
                   vtkInformationVector* outputVector);
   virtual int FillInputPortInformation(int port, vtkInformation* info);
   virtual int FillOutputPortInformation(int port, vtkInformation* info);

   /// Base ///
   int vtkDataSetToStdVector(vtkDataSet* input);
   int check(vtkDataSet* input);
   bool isCoincident(double p1[], double p2[]);

   /// ContourTree ///
   void getTree();
   void updateTree();
   NodeType getNodeType(int id);
   NodeType getNodeType(int id, TreeType type, MergeTree* tree);
   void getCriticalPoints();
   void clearTree();

   /// Skeleton ///
   void getSkeleton();
   void clearSkeleton();
   void getSkeletonNodes();
   void getSkeletonArcs();
   int getSkeletonScalars(const vector<double>& scalars,
                          vector<vector<double> >& skeletonScalars) const;

   /// Segmentation ///
   void getSegmentation(vtkDataSet* input);
   void clearSegmentation();

   /// Persistence Curve ///
   void getPersistenceCurve(TreeType type);
   void getCurves();

   /// Simplification Threshold Curve
   void getSimplificationThresholdCurve();

   /// Persistence Diagram ///
   void getPersistenceDiagramInfo(TreeType type, int pair, bool extremum, int& vertexId,
                                  int& nodeId, int& nodeType);
   void getPersistenceDiagram(TreeType type);
   void getDiagrams();

   int sample(unsigned int samplingLevel);

   int computeBarycenters();
   void computeSkeleton(unsigned int arcRes);
   void smoothSkeleton(unsigned int skeletonSmoothing);
   void smooth(const int idArc, bool order);

  private:
   /// Base ///
   bool UseAllCores;
   int ThreadNumber;
   int FieldId;
   bool isChecked_;
   bool isLoaded_;
   bool calculSegmentation_;
   bool lessPartition_;
   MergeTree* tree_;
   ParallelContourTree* contourTree_;
   vtkPolyData* skeletonNodes_;
   vtkPolyData* skeletonArcs_;
   vtkDataSet* segmentation_;
   vtkTable* CTPersistenceCurve_;
   vtkTable* MTPersistenceCurve_;
   vtkTable* STPersistenceCurve_;
   vtkTable* simplificationThresholdCurve_;
   vtkUnstructuredGrid* CTPersistenceDiagram_;
   vtkUnstructuredGrid* MTPersistenceDiagram_;
   vtkUnstructuredGrid* STPersistenceDiagram_;

   /// Void ///
   vtkUnstructuredGrid* voidUnstructuredGrid_;
   vtkPolyData* voidPolyData_;
   vtkTable* voidTable_;

   /// Configuration ///
   bool checkConformance_;
   string inputOffsetScalarFieldName_;
   bool useInputOffsetScalarField_;
   bool varyingMeshGeometry_;
   bool varyingMeshConnectivity_;
   bool varyingDataValues_;
   TreeType treeType_;
   string scalarField_;
   bool showMin_;
   bool showMax_;
   bool showSaddle1_;
   bool showSaddle2_;
   bool showArc_;
   unsigned int arcResolution_;
   int partitionNum_;
   unsigned int skeletonSmoothing_;
   int simplificationType_;
   double simplificationThreshold_;
   double simplificationThresholdBuffer_;
   bool simplifyPersistenceDiagram_;

   /// Computation handles
   bool toUpdateVertexSoSoffsets_;
   bool toComputeContourTree_;
   bool toComputeSimplification_;
   bool toUpdateTree_;
   bool toComputeSkeleton_;
   bool toComputeSegmentation_;
   bool toShowPersistenceCurve_;
   bool toShowPersistenceDiagram_;

   /// Convenient storage ///
   double deltaScalar_;
   unsigned int originalNumberOfCriticalPoints_;
   unsigned int numberOfVertices_;
   vtkTriangulation* triangulation_;
   vector<vector<int>>* vertexNeighbors_;
   vector<vector<double>>* vertexPositions_;
   vector<int>* vertexSoSoffsets_;
   vector<int>* criticalPoints_;
   vector<double>* vertexScalars_;
   vector<vector<double>>* inputScalars_;
   vector<string>* inputScalarsName_;
   vector<pair<pair<int, int>, double>>* mergePairs_;
   vector<pair<pair<int, int>, double>>* splitPairs_;
   vector<pair<pair<int, int>, double>>* pairs_;

   // treeType, SuperArc, several vertices list.
   vector<vector<vector<vector<int>>>>*    samples_;
   vector<vector<vector<vector<double>>>>* barycenters_;

   /// base code features ///
   bool needsToAbort();
   int updateProgress(const float& progress);
   int doIt(vtkDataSet* input, vtkPolyData* outputSkeletonNodes, vtkPolyData* outputSkeletonArcs,
            vtkDataSet* outputSegmentation, vtkTable* outputPersistenceCurve,
            vtkTable* outputSimplificationThresholdCurve,
            vtkUnstructuredGrid* outputPersistenceDiagram);
};

#endif  // _VTK_CONTOURTREE_H
