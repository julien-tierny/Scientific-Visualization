/* 
 * file:                  vtkFiberSurfaceFilter.cpp
 * description:           VTK wrapper for the FiberSurface package.
 * author(s):             Pavol Klacansky <pavol@klacansky.com>.
 *                        Julien Tierny <julien.tierny@lip6.fr>
 * date:                  March 2015.
 */

#include "extract_surface.hpp"
#include "extract_surface_grey.hpp"
#include "extract_surface_grey_minimal.hpp"

#include <vtkSetGet.h>

#include                  "vtkFiberSurfaceFilter.h"

typedef void (*extract_func)(const vtkIdType *tetIds,
                             const float *pointSet,
                             vtkPolyData *output,
                             const void *field0,
                             const void *field1,
                             const double *o,
                             const double *n,
                             const double &distance,
                             const double &d_length,
                             const double &cur_length,
                             const bool &visibleFibers,
                             vtkFloatArray *textureCoordinates,
                             vtkPoints *outputPoints,
                             const int &edgeId,
                             vtkIntArray *edgeIds);

/*
    resolve the function pointer once only to reduce code complexity
    we can swap implementations on run time (simple, grey, and grey minimal)
    performance seems same (even it is not inlined)

    TODO Pavol does not know how the default arguments are handled in assembly
    so there might be a bug

	HELL IS REAL
*/
static extract_func
get_extract_func(const int dataType0, const int dataType1, const int algorithm)
{
  extract_func fp;
  switch(dataType0){
        case VTK_CHAR:
          switch(dataType1){
            vtkTemplateMacro(({
              if (algorithm == vtkFiberSurfaceFilter::SIMPLE)
                fp = extract_surface<char, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY)
			    fp = extract_surface_grey<char, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY_MINIMAL)
                fp = extract_surface_grey_minimal<char, VTK_TT>;
              else
                std::cerr << "[vtkFiberSurfaceFilter] Unsupported algorithm type\n";
              })
            );
          }
          break;
        case VTK_DOUBLE:
          switch(dataType1){
            vtkTemplateMacro(({
              if (algorithm == vtkFiberSurfaceFilter::SIMPLE)
                fp = extract_surface<double, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY)
                fp = extract_surface_grey<double, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY_MINIMAL)
                fp = extract_surface_grey_minimal<double, VTK_TT>;
              else
                std::cerr << "[vtkFiberSurfaceFilter] Unsupported algorithm type\n";
              })
            );
          }
          break;
        case VTK_FLOAT:
          switch(dataType1){
            vtkTemplateMacro(({
              if (algorithm == vtkFiberSurfaceFilter::SIMPLE)
                fp = extract_surface<float, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY)
                fp = extract_surface_grey<float, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY_MINIMAL)
                fp = extract_surface_grey_minimal<float, VTK_TT>;
              else
                std::cerr << "[vtkFiberSurfaceFilter] Unsupported algorithm type\n";
              })
            );
          } 
          break;
        case VTK_INT:
          switch(dataType1){
            vtkTemplateMacro(({
              if (algorithm == vtkFiberSurfaceFilter::SIMPLE)
                fp = extract_surface<int, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY)
                fp = extract_surface_grey<int, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY_MINIMAL)
                fp = extract_surface_grey_minimal<int, VTK_TT>;
              else
                std::cerr << "[vtkFiberSurfaceFilter] Unsupported algorithm type\n";
              })
            );
          }
          break;
        case VTK_UNSIGNED_CHAR:
          switch(dataType1){
            vtkTemplateMacro(({
              if (algorithm == vtkFiberSurfaceFilter::SIMPLE)
                fp = extract_surface<unsigned char, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY)
                fp = extract_surface_grey<unsigned char, VTK_TT>;
              else if (algorithm == vtkFiberSurfaceFilter::GREY_MINIMAL)
                fp = extract_surface_grey_minimal<unsigned char, VTK_TT>;
              else
                std::cerr << "[vtkFiberSurfaceFilter] Unsupported algorithm type\n";
              })
            );
          }
          break;
        default:
          {
            std::cerr << "[vtkFiberSurfaceFilter] Unsupported data type :(\n";
          }
          break;
      }
 return fp;
}


vtkStandardNewMacro(vtkFiberSurfaceFilter)

vtkFiberSurfaceFilter::vtkFiberSurfaceFilter(){
  
  input = nullptr;
  output = nullptr;

  SetNumberOfInputPorts(2);
  
  ThreadBalance = 1;
  ThreadStrategy = DOMAIN_TETRAHEDRA;
  ThreadNumber = threadNumber_;
  Manifold = false;
  Implementation = REGULAR;
  Algorithm = SIMPLE;

  BVH_ = NULL;
  octree_ = NULL;
  ShowCandidates = false;
  VisibleFibers = false;
  PolygonSegmentation = false;
  BvhMinimumTetNumber = 8;
}

vtkFiberSurfaceFilter::~vtkFiberSurfaceFilter(){

  if(BVH_){
    bvh_free(BVH_);
  }
  
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
}

int vtkFiberSurfaceFilter::FillInputPortInformation(int port, 
  vtkInformation *info){

  switch (port) {
    case 0:
      info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),
                "vtkUnstructuredGrid");
      return 1;
    case 1:
      info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),
                "vtkUnstructuredGrid");
      return 1;
  }

  return 0;  
}

template<class dataType0, class dataType1> int vtkFiberSurfaceFilter::buildBVH(
  vtkUnstructuredGrid *input, void *pU, void *pV){

  Timer t;
  
  uint32_t nodes_per_leaf = BvhMinimumTetNumber;
  
  if(nodes_per_leaf < 1)
    nodes_per_leaf = 1;
  if(nodes_per_leaf > 16)
    nodes_per_leaf = 16;
  
  int cellNumber = input->GetNumberOfCells();
  
  vtkIdType *inputCells = input->GetCells()->GetPointer();
  
  float *bboxs = new float[4*cellNumber];
  
  dataType0 *u = (dataType0 *) pU;
  dataType1 *v = (dataType1 *) pV;
  
  // init -- openMP... slight gain 
#ifdef withOpenMP
#pragma omp parallel for num_threads(ThreadNumber)
#endif
  for(int i = 0; i < cellNumber; i++){
    
    float *bbox = bboxs + 4*i;
    bbox[0] = bbox[1] = FLT_MAX;
    bbox[2] = bbox[3] = -FLT_MAX;
  
    vtkIdType *cell = &(inputCells[5*i + 1]);
    
    for(int j = 0; j < 4; j++){
      bbox[0] = fminf(bbox[0], u[cell[j]]);
      bbox[1] = fminf(bbox[1], v[cell[j]]);
      bbox[2] = fmaxf(bbox[2], u[cell[j]]);
      bbox[3] = fmaxf(bbox[3], v[cell[j]]);

    }
  }
  
  BVH_ = bvh_alloc(cellNumber, nodes_per_leaf);
  
  bvh_build(BVH_, bboxs);
  bvh_stats(BVH_);
  {
    stringstream msg;
    msg << "[vtkFiberSurfaceFilter] BVH pre-process time: "
      << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), 2);
  }
  
  delete[] bboxs;
  
  return 0;
}

template<class dataType0, class dataType1> 
  int vtkFiberSurfaceFilter::buildOctree(
    vtkUnstructuredGrid *input, void *pU, void *pV){
  
  RangeDrivenOctree<dataType0, dataType1> *octree = 
    new RangeDrivenOctree<dataType0, dataType1>;
    
  octree->setWrapper(this);
  octree->setVertexNumber(input->GetNumberOfPoints());
  octree->setCellNumber(input->GetNumberOfCells());
  octree->setCellList(input->GetCells()->GetPointer());
  octree->setPointList((float *) input->GetPoints()->GetVoidPointer(0));
  octree->setRange((dataType0 *) pU, (dataType1 *) pV);
  octree->setLeafMinimumRangeAreaRatio(octreeMinimumRangeAreaRatio_);
  
  octree->build();
  
  octree_ = (void *) octree;
    
  return 0;
}

double vtkFiberSurfaceFilter::doIt(vtkUnstructuredGrid *input,
  vtkUnstructuredGrid *polygon, vtkPolyData *output){

  threadedOutput_.clear();
  threadedPoints_.clear();
  threadedCells_.clear();
  threadedTextureCoordinates_.clear();
  threadedEdgeIds_.clear();
  
  // allocate output

  appender_ = vtkSmartPointer<vtkAppendPolyData>::New();
  threadedOutput_.resize(ThreadNumber);
  threadedPoints_.resize(ThreadNumber);
  threadedCells_.resize(ThreadNumber);
  threadedTextureCoordinates_.resize(ThreadNumber);
  threadedEdgeIds_.resize(ThreadNumber);
  
  for(int i = 0; i < (int) threadedOutput_.size(); i++){
    threadedOutput_[i] = vtkSmartPointer<vtkPolyData>::New();
    threadedPoints_[i] = vtkSmartPointer<vtkPoints>::New();
    threadedCells_[i] = vtkSmartPointer<vtkCellArray>::New();
    if(VisibleFibers)
      threadedTextureCoordinates_[i] = vtkSmartPointer<vtkFloatArray>::New();
    if(PolygonSegmentation){
      threadedEdgeIds_[i] = vtkSmartPointer<vtkIntArray>::New();
      threadedEdgeIds_[i]->SetNumberOfComponents(1);
      threadedEdgeIds_[i]->SetName("Polygon Edge Id");
      threadedOutput_[i]->GetCellData()->AddArray(threadedEdgeIds_[i]);
    }
    threadedOutput_[i]->SetPoints(threadedPoints_[i]);
    threadedOutput_[i]->SetPolys(threadedCells_[i]);
    
    if(VisibleFibers){
      threadedTextureCoordinates_[i]->SetNumberOfComponents(2);
      threadedOutput_[i]->GetPointData()->SetTCoords(
        threadedTextureCoordinates_[i]);
    }
  }
  
  // retrieve info about the input polygon  
  vtkPoints *line_points = polygon->GetPoints();

  vtkDataArray *polygonUComponent = NULL;
  vtkDataArray *polygonVComponent = NULL;
  if(PolygonUComponent.size()){
    polygonUComponent = polygon->GetPointData()->GetArray(
      PolygonUComponent.data());
  }
  if(PolygonVComponent.size()){
    polygonVComponent = polygon->GetPointData()->GetArray(
      PolygonVComponent.data());
  }
  
  // retrieve info about the input data-set
  vtkDataArray *fields[2] = { NULL, NULL };
  
  if(DataUComponent.size()){
    fields[0] = input->GetPointData()->GetArray(DataUComponent.data());
  }
  if(DataVComponent.size()){
    fields[1] = input->GetPointData()->GetArray(DataVComponent.data());
  }
  
  // back up option
  if(!fields[0]){
    fields[0] = input->GetPointData()->GetArray(0);
  }
  if(!fields[1]){
    fields[1] = input->GetPointData()->GetArray(1);
  }

  double p[3];
  double d[2], n[2], o[2];

  const vtkIdType polygon_length = polygon->GetNumberOfCells();

  // hold the length of the polyline
  // TODO makes sense if and only if the input is a polyline
  double total_length = 0.0;
  double cur_length = 0.0;

  int dataType0 = fields[0]->GetDataType();
  int dataType1 = fields[1]->GetDataType();
  extract_func fp = get_extract_func(dataType0, dataType1, Algorithm);
  
  void *field0 = fields[0]->GetVoidPointer(0);
  void *field1 = fields[1]->GetVoidPointer(0);
  
  float *pointSet = (float *) input->GetPoints()->GetVoidPointer(0);
  
  for (vtkIdType i = 0; i != polygon_length; ++i) {
    
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Processing polygon edge #" << i << endl;
      dMsg(cout, msg.str(), 3);
    }

    vtkCell *line = polygon->GetCell(i);

    line_points->GetPoint(line->GetPointId(0), p);
    if(polygonUComponent){
      polygonUComponent->GetTuple(line->GetPointId(0), &(p[0]));
    }
    if(polygonVComponent){
      polygonVComponent->GetTuple(line->GetPointId(0), &(p[1]));
    }
    o[0] = p[0];
    o[1] = p[1];

    line_points->GetPoint(line->GetPointId(1), p);
    if(polygonUComponent){
      polygonUComponent->GetTuple(line->GetPointId(1), &(p[0]));
    }
    if(polygonVComponent){
      polygonVComponent->GetTuple(line->GetPointId(1), &(p[1]));
    }   
    d[0] = p[0] - o[0];
    d[1] = p[1] - o[1];

    
    // test if the line is not of zero length
    if (d[0]*d[0] + d[1]*d[1] <= FLT_EPSILON) {
      std::cerr << name << ": Zero length line segment; skipping\n";
      continue;
    }

    n[0] = d[1];
    n[1] = -d[0];
    const double d_length = std::sqrt(n[0]*n[0] + n[1]*n[1]);
    total_length += d_length;
    n[0] /= d_length;
    n[1] /= d_length;

    // compute line segment normal and line equation
    const double distance = n[0]*o[0] + n[1]*o[1];

    // extract fiber surface
    const vtkIdType cells_length = input->GetNumberOfCells();
    
    vtkIdType *inputCells = input->GetCells()->GetPointer();
#ifdef withOpenMP
#pragma omp parallel for num_threads(ThreadNumber)
#endif
    for(int j = 0; j < cells_length; j++){
    
      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
#endif
      
      fp(&inputCells[5*j + 1],
        pointSet, threadedOutput_[threadId].GetPointer(), 
        field0, field1, 
        o, n, distance, d_length, cur_length, VisibleFibers,
        threadedTextureCoordinates_[threadId].GetPointer(),
        threadedPoints_[threadId].GetPointer(), 
        i, threadedEdgeIds_[threadId]);
        
    }
    
    cur_length = total_length;
  }
  
  reduce();
  
  return total_length;
}

double vtkFiberSurfaceFilter::doItOctree(vtkUnstructuredGrid *input, 
  vtkUnstructuredGrid *polygon, vtkPolyData *output){

  int polygonEdgeThreadNumber = 1;
  int fiberSurfaceThreadNumber = ThreadNumber;
  
  threadedOutput_.clear();
  threadedPoints_.clear();
  threadedCells_.clear();
  threadedTextureCoordinates_.clear();
  threadedEdgeIds_.clear();
  
  // allocate output

  appender_ = vtkSmartPointer<vtkAppendPolyData>::New();
  threadedOutput_.resize(ThreadNumber);
  threadedPoints_.resize(ThreadNumber);
  threadedCells_.resize(ThreadNumber);
  threadedTextureCoordinates_.resize(ThreadNumber);
  threadedEdgeIds_.resize(ThreadNumber);
  
  for(int i = 0; i < (int) threadedOutput_.size(); i++){
    threadedOutput_[i] = vtkSmartPointer<vtkPolyData>::New();
    threadedPoints_[i] = vtkSmartPointer<vtkPoints>::New();
    threadedCells_[i] = vtkSmartPointer<vtkCellArray>::New();
    if(PolygonSegmentation){
      threadedEdgeIds_[i] = vtkSmartPointer<vtkIntArray>::New();
      threadedEdgeIds_[i]->SetNumberOfComponents(1);
      threadedEdgeIds_[i]->SetName("Polygon Edge Id");
      threadedOutput_[i]->GetCellData()->AddArray(threadedEdgeIds_[i]);
    }
    if(VisibleFibers)
      threadedTextureCoordinates_[i] = vtkSmartPointer<vtkFloatArray>::New();
    threadedOutput_[i]->SetPoints(threadedPoints_[i]);
    threadedOutput_[i]->SetPolys(threadedCells_[i]);
    
    if(VisibleFibers){
      threadedTextureCoordinates_[i]->SetNumberOfComponents(2);
      threadedOutput_[i]->GetPointData()->SetTCoords(
        threadedTextureCoordinates_[i]);
    }
  }
  
  // retrieve info about the input polygon  
  vtkPoints *line_points = polygon->GetPoints();

  vtkDataArray *polygonUComponent = NULL;
  vtkDataArray *polygonVComponent = NULL;
  if(PolygonUComponent.size()){
    polygonUComponent = polygon->GetPointData()->GetArray(
      PolygonUComponent.data());
  }
  if(PolygonVComponent.size()){
    polygonVComponent = polygon->GetPointData()->GetArray(
      PolygonVComponent.data());
  }
  
  // retrieve info about the input data-set
  vtkDataArray *fields[2] = { NULL, NULL };
  
  if(DataUComponent.size()){
    fields[0] = input->GetPointData()->GetArray(DataUComponent.data());
  }
  if(DataVComponent.size()){
    fields[1] = input->GetPointData()->GetArray(DataVComponent.data());
  }
  
  // back up option
  if(!fields[0]){
    fields[0] = input->GetPointData()->GetArray(0);
  }
  if(!fields[1]){
    fields[1] = input->GetPointData()->GetArray(1);
  }



  const vtkIdType polygon_length = polygon->GetNumberOfCells();

  // hold the length of the polyline
  double total_length = 0.0;
  double cur_length = 0.0;

  int dataType0 = fields[0]->GetDataType();
  int dataType1 = fields[1]->GetDataType();
  
  extract_func fp = get_extract_func(dataType0, dataType1, Algorithm);
  
  void *field0 = fields[0]->GetVoidPointer(0);
  void *field1 = fields[1]->GetVoidPointer(0);
  
  float *pointSet = (float *) input->GetPoints()->GetVoidPointer(0);
  
  // extract fiber surface
  vtkIdType *inputCells = input->GetCells()->GetPointer();
  
  switch(ThreadStrategy){
    case DOMAIN_TETRAHEDRA:
      polygonEdgeThreadNumber = 1;
      fiberSurfaceThreadNumber = ThreadNumber;
      break;
    case POLYGON_EDGES:
      polygonEdgeThreadNumber = ThreadNumber;
      fiberSurfaceThreadNumber = 1;
      break;
    case AUTOMATIC:
      // if you have more edges than threads, let's go full polygon edges
      if(polygon_length > ThreadBalance*ThreadNumber){
              polygonEdgeThreadNumber = ThreadNumber;
      fiberSurfaceThreadNumber = 1;
      }
      break;
  }

  {
    stringstream msg;
    msg << "[vtkFiberSurfaceFilter] Using "
      << polygonEdgeThreadNumber << " threads for polygon edges..." 
      << endl;
    msg << "[vtkFiberSurfaceFilter] Using "
      << fiberSurfaceThreadNumber
      << " threads for fiber surfaces..."
      << endl;
    dMsg(cout, msg.str(), 2);
  }
 
  // octree query data
  vector<vector<int> > cellIds(ThreadNumber);
  
  vector<vector<int> > candidateTets(polygon_length);
 
#ifdef withOpenMP
#pragma omp parallel for num_threads(polygonEdgeThreadNumber)
#endif
  for(int i = 0; i < polygon_length; i++){
   
    int edgeThreadId = 0;
#ifdef withOpenMP
    edgeThreadId = omp_get_thread_num();
#endif
    
    double p[3];
    double d[2], n[2], o[2];
    
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Processing polygon edge #" << i 
        << "..." << endl;
      dMsg(cout, msg.str(), 3);
    }

    vtkCell *line = polygon->GetCell(i);

    line_points->GetPoint(line->GetPointId(0), p);
    if(polygonUComponent){
      polygonUComponent->GetTuple(line->GetPointId(0), &(p[0]));
    }
    if(polygonVComponent){
      polygonVComponent->GetTuple(line->GetPointId(0), &(p[1]));
    }
    o[0] = p[0];
    o[1] = p[1];

    line_points->GetPoint(line->GetPointId(1), p);
    if(polygonUComponent){
      polygonUComponent->GetTuple(line->GetPointId(1), &(p[0]));
    }
    if(polygonVComponent){
      polygonVComponent->GetTuple(line->GetPointId(1), &(p[1]));
    }   
    d[0] = p[0] - o[0];
    d[1] = p[1] - o[1];

    
    // test if the line is not of zero length
    if (d[0]*d[0] + d[1]*d[1] <= FLT_EPSILON) {
      std::cerr << name << ": Zero length line segment; skipping\n";
      continue;
    }

    n[0] = d[1];
    n[1] = -d[0];
    const double d_length = std::sqrt(n[0]*n[0] + n[1]*n[1]);
    total_length += d_length;
    n[0] /= d_length;
    n[1] /= d_length;

    // compute line segment normal and line equation
    const double distance = n[0]*o[0] + n[1]*o[1];

    // query the octree
    Timer query;
    
    switch(dataType0){
      case VTK_CHAR:
        switch(dataType1){
          vtkTemplateMacro(
            (
            ((RangeDrivenOctree<char, VTK_TT> *) octree_)->rangeSegmentQuery(
              o, p, cellIds[edgeThreadId])
            )
          );
        }
        break;
      case VTK_DOUBLE:
        switch(dataType1){
          vtkTemplateMacro(
            (
            ((RangeDrivenOctree<double, VTK_TT> *) octree_)->rangeSegmentQuery(
              o, p, cellIds[edgeThreadId])
            )
          );
        }
        break;
      case VTK_FLOAT:
        switch(dataType1){
          vtkTemplateMacro(
            (
            ((RangeDrivenOctree<float, VTK_TT> *) octree_)->rangeSegmentQuery(
              o, p, cellIds[edgeThreadId])
            )
          );
        }
        break;
      case VTK_INT:
        switch(dataType1){
          vtkTemplateMacro(
            (
            ((RangeDrivenOctree<int, VTK_TT> *) octree_)->rangeSegmentQuery(
              o, p, cellIds[edgeThreadId])
            )
          );
        }
        break;
      case VTK_UNSIGNED_CHAR:
        switch(dataType1){
          vtkTemplateMacro(
            (
            ((RangeDrivenOctree<unsigned char, VTK_TT> *) 
              octree_)->rangeSegmentQuery(o, p, cellIds[edgeThreadId])
            )
          );
        }
        break;
      default:
        {
          stringstream msg;
          msg << "[vtkFiberSurfaceFilter] Unsupported data type :(" << endl;
          dMsg(cerr, msg.str(), 1);
        }
        break;
    }
    
    if(ShowCandidates){
      candidateTets[i] = cellIds[edgeThreadId];
    }
    
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Octree query (edge #"
        << i << ") done in "
        << query.getElapsedTime() << " s. (" 
        << cellIds[edgeThreadId].size()
        << " tets)"<< endl;
      dMsg(cout, msg.str(), 2);
    }
    
    Timer extraction;
#ifdef withOpenMP
#pragma omp parallel for num_threads(fiberSurfaceThreadNumber)
#endif
    for(int j = 0; j < (int) cellIds[edgeThreadId].size(); j++){
    
      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
      if(fiberSurfaceThreadNumber == 1)
        threadId = edgeThreadId;
#endif
      
      fp( &inputCells[5*cellIds[edgeThreadId][j] + 1], 
        pointSet, threadedOutput_[threadId].GetPointer(), 
        field0, field1, 
        o, n, distance, d_length, cur_length, VisibleFibers,
        threadedTextureCoordinates_[threadId].GetPointer(),
        threadedPoints_[threadId].GetPointer(),
        i, threadedEdgeIds_[threadId]);
      
    }
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Polygon-edge #"
        << i << " fiber-processed in: "
        << extraction.getElapsedTime() << " s." << endl;
      dMsg(cout, msg.str(), 2);
    }
    
    cur_length = total_length;
  }
  
  reduce();
  
  if(ShowCandidates){
    
    Timer candidateTimer;
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Preparing candidate tet list..."
        << endl;
      dMsg(cout, msg.str(), 4);
    }
   
    int nonCandidateCellNumber = output->GetNumberOfCells();
   
    vtkSmartPointer<vtkIntArray> edgeCandidate =   
      vtkSmartPointer<vtkIntArray>::New();
    edgeCandidate->SetNumberOfComponents(1);
    edgeCandidate->SetName("CandidateEdgeId");
    edgeCandidate->SetNumberOfTuples(nonCandidateCellNumber);
    output->GetCellData()->AddArray(edgeCandidate);

    int edgeId = -1;
    long long int *cells = input->GetCells()->GetPointer();
    for(int i = 0; i < edgeCandidate->GetNumberOfTuples(); i++){
      edgeCandidate->SetTuple1(i, edgeId);
    }
    
    for(int i = 0; i < (int) candidateTets.size(); i++){
      edgeId = i;
      
      for(int j = 0; j < (int) candidateTets[i].size(); j++){
        
        long long int *cell = &(cells[5*candidateTets[i][j] + 1]);
        
        double p[3];
        vector<int> newIds(4);
        for(int k = 0; k < 4; k++){
          input->GetPoint(cell[k], p);
          newIds[k] = output->GetPoints()->InsertNextPoint(p);
        }
        
        vtkSmartPointer<vtkIdList> vertexList = 
          vtkSmartPointer<vtkIdList>::New();
        vertexList->SetNumberOfIds(3);
        for(int k = 0; k < 3; k++){
          for(int l = 0; l < 3; l++){
            vertexList->SetId(l, newIds[(k + l)%4]);
          }
          edgeCandidate->InsertNextTuple1(edgeId);
          if(PolygonSegmentation){
            output->GetCellData()->GetArray(
              "Polygon Edge Id")->InsertNextTuple1(edgeId);
          }
          output->GetPolys()->InsertNextCell(vertexList);
        }
        
        // last triangle
        vertexList->SetId(0, newIds[0]);
        vertexList->SetId(1, newIds[1]);
        vertexList->SetId(2, newIds[3]);
        edgeCandidate->InsertNextTuple1(edgeId);
        if(PolygonSegmentation){
          output->GetCellData()->GetArray(
              "Polygon Edge Id")->InsertNextTuple1(edgeId);
        }
        output->GetPolys()->InsertNextCell(vertexList);
      }
    }
    
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Candidate tets prepared in "
        << candidateTimer.getElapsedTime() << " s." << endl;
      dMsg(cout, msg.str(), 2);
    }
  }
  
  return total_length;
}

double vtkFiberSurfaceFilter::doItBVH(vtkUnstructuredGrid *input, 
  vtkUnstructuredGrid *polygon, vtkPolyData *output){

  
  int polygonEdgeThreadNumber = 1;
  int fiberSurfaceThreadNumber = ThreadNumber;
 

  
  threadedOutput_.clear();
  threadedPoints_.clear();
  threadedCells_.clear();
  threadedTextureCoordinates_.clear();
  threadedEdgeIds_.clear();
  threadedEdgeIds_.resize(ThreadNumber);
  
  // allocate output

  appender_ = vtkSmartPointer<vtkAppendPolyData>::New();
  threadedOutput_.resize(ThreadNumber);
  threadedPoints_.resize(ThreadNumber);
  threadedCells_.resize(ThreadNumber);
  threadedTextureCoordinates_.resize(ThreadNumber);
  
  for(int i = 0; i < (int) threadedOutput_.size(); i++){
    threadedOutput_[i] = vtkSmartPointer<vtkPolyData>::New();
    threadedPoints_[i] = vtkSmartPointer<vtkPoints>::New();
    threadedCells_[i] = vtkSmartPointer<vtkCellArray>::New();
    if(PolygonSegmentation){
      threadedEdgeIds_[i] = vtkSmartPointer<vtkIntArray>::New();
      threadedEdgeIds_[i]->SetNumberOfComponents(1);
      threadedEdgeIds_[i]->SetName("Polygon Edge Id");
      threadedOutput_[i]->GetCellData()->AddArray(threadedEdgeIds_[i]);
    }
    if(VisibleFibers)
      threadedTextureCoordinates_[i] = vtkSmartPointer<vtkFloatArray>::New();
    threadedOutput_[i]->SetPoints(threadedPoints_[i]);
    threadedOutput_[i]->SetPolys(threadedCells_[i]);
    
    if(VisibleFibers){
      threadedTextureCoordinates_[i]->SetNumberOfComponents(2);
      threadedOutput_[i]->GetPointData()->SetTCoords(
        threadedTextureCoordinates_[i]);
    }
  }
  
  // retrieve info about the input polygon  
  vtkPoints *line_points = polygon->GetPoints();

  vtkDataArray *polygonUComponent = NULL;
  vtkDataArray *polygonVComponent = NULL;
  if(PolygonUComponent.size()){
    polygonUComponent = polygon->GetPointData()->GetArray(
      PolygonUComponent.data());
  }
  if(PolygonVComponent.size()){
    polygonVComponent = polygon->GetPointData()->GetArray(
      PolygonVComponent.data());
  }
  
  // retrieve info about the input data-set
  vtkDataArray *fields[2] = { NULL, NULL };
  
  if(DataUComponent.size()){
    fields[0] = input->GetPointData()->GetArray(DataUComponent.data());
  }
  if(DataVComponent.size()){
    fields[1] = input->GetPointData()->GetArray(DataVComponent.data());
  }
  
  // back up option
  if(!fields[0]){
    fields[0] = input->GetPointData()->GetArray(0);
  }
  if(!fields[1]){
    fields[1] = input->GetPointData()->GetArray(1);
  }



  const vtkIdType polygon_length = polygon->GetNumberOfCells();

  // hold the length of the polyline
  double total_length = 0.0;
  double cur_length = 0.0;

  int dataType0 = fields[0]->GetDataType();
  int dataType1 = fields[1]->GetDataType();
 
  extract_func fp = get_extract_func(dataType0, dataType1, Algorithm);
  
  void *field0 = fields[0]->GetVoidPointer(0);
  void *field1 = fields[1]->GetVoidPointer(0);
  
  float *pointSet = (float *) input->GetPoints()->GetVoidPointer(0);
  
  // extract fiber surface
  vtkIdType *inputCells = input->GetCells()->GetPointer();

  switch(ThreadStrategy){
    case DOMAIN_TETRAHEDRA:
      polygonEdgeThreadNumber = 1;
      fiberSurfaceThreadNumber = ThreadNumber;
      break;
    case POLYGON_EDGES:
      polygonEdgeThreadNumber = ThreadNumber;
      fiberSurfaceThreadNumber = 1;
      break;
    case AUTOMATIC:
      // if you have more edges than threads, let's go full polygon edges
      if(polygon_length > ThreadBalance*ThreadNumber){
              polygonEdgeThreadNumber = ThreadNumber;
      fiberSurfaceThreadNumber = 1;
      }
      break;
  }

  {
    stringstream msg;
    msg << "[vtkFiberSurfaceFilter] Using "
      << polygonEdgeThreadNumber << " threads for polygon edges..." 
      << endl;
    msg << "[vtkFiberSurfaceFilter] Using "
      << fiberSurfaceThreadNumber
      << " threads for fiber surfaces..."
      << endl;
    dMsg(cout, msg.str(), 2);
  }
 
  // bvh query data
  const uint32_t ids_capacity = 64*1024*1024;
  vector<uint32_t *> threadedIdList(ThreadNumber);
  for(int i = 0; i < (int) threadedIdList.size(); i++){
    threadedIdList[i] = new uint32_t[ids_capacity];
  }
  
  vector<vector<int> > candidateTets(polygon_length);
 
#ifdef withOpenMP
#pragma omp parallel for num_threads(polygonEdgeThreadNumber)
#endif
  for(int i = 0; i < polygon_length; i++) {
    
    int edgeThreadId = 0;
#ifdef withOpenMP
    edgeThreadId = omp_get_thread_num();
#endif
    
    double p[3];
    double d[2], n[2], o[2];
    
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Processing polygon edge #" << i 
        << "..." << endl;
      dMsg(cout, msg.str(), 3);
    }

    vtkCell *line = polygon->GetCell(i);

    line_points->GetPoint(line->GetPointId(0), p);
    if(polygonUComponent){
      polygonUComponent->GetTuple(line->GetPointId(0), &(p[0]));
    }
    if(polygonVComponent){
      polygonVComponent->GetTuple(line->GetPointId(0), &(p[1]));
    }
    o[0] = p[0];
    o[1] = p[1];

    line_points->GetPoint(line->GetPointId(1), p);
    if(polygonUComponent){
      polygonUComponent->GetTuple(line->GetPointId(1), &(p[0]));
    }
    if(polygonVComponent){
      polygonVComponent->GetTuple(line->GetPointId(1), &(p[1]));
    }   
    d[0] = p[0] - o[0];
    d[1] = p[1] - o[1];

    
    // test if the line is not of zero length
    if (d[0]*d[0] + d[1]*d[1] <= FLT_EPSILON) {
      std::cerr << name << ": Zero length line segment; skipping\n";
      continue;
    }

    n[0] = d[1];
    n[1] = -d[0];
    const double d_length = std::sqrt(n[0]*n[0] + n[1]*n[1]);
    total_length += d_length;
    n[0] /= d_length;
    n[1] /= d_length;

    // compute line segment normal and line equation
    const double distance = n[0]*o[0] + n[1]*o[1];

    // query the BVH
    float l[4];
    l[0] = o[0];
    l[1] = o[1];
    l[2] = p[0];
    l[3] = p[1];
   
    Timer query;
    const uint32_t ids_length = bvh_intersect(BVH_, 
      threadedIdList[edgeThreadId], ids_capacity, l);
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] BVH query (edge #"
        << i << ") done in "
        << query.getElapsedTime() << " s. (" << (int) ids_length
        << " tets)"<< endl;
      dMsg(cout, msg.str(), 2);
    }
    
    Timer extraction;
#ifdef withOpenMP
#pragma omp parallel for num_threads(fiberSurfaceThreadNumber)
#endif
    for(int j = 0; j < (int) ids_length; j++){
    
      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
      if(fiberSurfaceThreadNumber == 1)
        threadId = edgeThreadId;
#endif
     
      fp(&inputCells[5*threadedIdList[edgeThreadId][j] + 1],
        pointSet, threadedOutput_[threadId].GetPointer(), 
        field0, field1, 
        o, n, distance, d_length, cur_length, VisibleFibers,
        threadedTextureCoordinates_[threadId].GetPointer(),
        threadedPoints_[threadId].GetPointer(),
        i, threadedEdgeIds_[threadId]);
      
    }
    
    if(ShowCandidates){
      candidateTets[i].resize(ids_length);
      for(int j = 0; j < (int) candidateTets[i].size(); j++){
        candidateTets[i][j] = threadedIdList[edgeThreadId][j];
      }
    }
    
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Polygon-edge #"
        << i << " fiber-processed in: "
        << extraction.getElapsedTime() << " s." << endl;
      dMsg(cout, msg.str(), 2);
    }
    
    cur_length = total_length;
  }
  
  reduce();
  
  if(ShowCandidates){
    
    Timer candidateTimer;
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Preparing candidate tet list..."
        << endl;
      dMsg(cout, msg.str(), 4);
    }
   
    int nonCandidateCellNumber = output->GetNumberOfCells();
   
    vtkSmartPointer<vtkIntArray> edgeCandidate =   
      vtkSmartPointer<vtkIntArray>::New();
    edgeCandidate->SetNumberOfComponents(1);
    edgeCandidate->SetName("CandidateEdgeId");
    edgeCandidate->SetNumberOfTuples(nonCandidateCellNumber);
    output->GetCellData()->AddArray(edgeCandidate);

    int edgeId = -1;
    long long int *cells = input->GetCells()->GetPointer();
    for(int i = 0; i < edgeCandidate->GetNumberOfTuples(); i++){
      edgeCandidate->SetTuple1(i, edgeId);
    }
    
    for(int i = 0; i < (int) candidateTets.size(); i++){
      edgeId = i;
      
      for(int j = 0; j < (int) candidateTets[i].size(); j++){
        
        long long int *cell = &(cells[5*candidateTets[i][j] + 1]);
        
        double p[3];
        vector<int> newIds(4);
        for(int k = 0; k < 4; k++){
          input->GetPoint(cell[k], p);
          newIds[k] = output->GetPoints()->InsertNextPoint(p);
        }
        
        vtkSmartPointer<vtkIdList> vertexList = 
          vtkSmartPointer<vtkIdList>::New();
        vertexList->SetNumberOfIds(3);
        for(int k = 0; k < 3; k++){
          for(int l = 0; l < 3; l++){
            vertexList->SetId(l, newIds[(k + l)%4]);
          }
          edgeCandidate->InsertNextTuple1(edgeId);
          if(PolygonSegmentation){
            output->GetCellData()->GetArray(
              "Polygon Edge Id")->InsertNextTuple1(edgeId);
          }
          output->GetPolys()->InsertNextCell(vertexList);
        }
        
        // last triangle
        vertexList->SetId(0, newIds[0]);
        vertexList->SetId(1, newIds[1]);
        vertexList->SetId(2, newIds[3]);
        edgeCandidate->InsertNextTuple1(edgeId);
        if(PolygonSegmentation){
          output->GetCellData()->GetArray(
              "Polygon Edge Id")->InsertNextTuple1(edgeId);
        }
        output->GetPolys()->InsertNextCell(vertexList);
      }
    }
    
    {
      stringstream msg;
      msg << "[vtkFiberSurfaceFilter] Candidate tets prepared in "
        << candidateTimer.getElapsedTime() << " s." << endl;
      dMsg(cout, msg.str(), 2);
    }
  }
  
  for(int i = 0; i < (int) threadedIdList.size(); i++){
    delete[] threadedIdList[i];
  }
  
  return total_length;
}

int vtkFiberSurfaceFilter::reduce() {

  Timer timer;
 
  vtkSmartPointer<vtkPoints> reducedPoints = 
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> reducedCells =
    vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyData> reducedSurface = 
    vtkSmartPointer<vtkPolyData>::New();

  int triangleNumber = 0;
  vector<int> triangleOffsets(ThreadNumber);
  
  for(int i = 0; i < ThreadNumber; i++){
    triangleOffsets[i] = triangleNumber;
    triangleNumber += threadedOutput_[i]->GetNumberOfCells();
  }
  
  reducedPoints->SetNumberOfPoints(3*triangleNumber);
  reducedCells->SetNumberOfCells(triangleNumber);
  
  // init surface
  reducedSurface->SetPoints(reducedPoints);
  reducedSurface->SetPolys(reducedCells);
  
  float *points = (float *) reducedPoints->GetVoidPointer(0);
  vtkSmartPointer<vtkIdTypeArray> cellList = 
    vtkSmartPointer<vtkIdTypeArray>::New();
  cellList->SetNumberOfTuples(4*triangleNumber);
  long long int *cells = (long long int *) cellList->GetVoidPointer(0);
 
  vtkSmartPointer<vtkFloatArray> reducedFiberCoordinates =
    vtkSmartPointer<vtkFloatArray>::New();
  if(VisibleFibers){
    reducedFiberCoordinates->SetNumberOfComponents(2);
    reducedFiberCoordinates->SetNumberOfTuples(3*triangleNumber);
  }
  float *reducedFiberPointer = 
    (float *) reducedFiberCoordinates->GetVoidPointer(0);
  
  vtkSmartPointer<vtkIntArray> reducedEdgeIds =
    vtkSmartPointer<vtkIntArray>::New();
  if(PolygonSegmentation){
    reducedEdgeIds->SetNumberOfComponents(1);
    reducedEdgeIds->SetNumberOfTuples(triangleNumber);
  }
  int *edgeIdsPointer = (int *) reducedEdgeIds->GetVoidPointer(0); 
    
#ifdef withOpenMP
#pragma omp parallel for num_threads(ThreadNumber)
#endif
  for(int i = 0; i < ThreadNumber; i++){
    
    float *threadedPoints = (float *) threadedPoints_[i]->GetVoidPointer(0);
    long long int *threadedCells = (long long int *)
      threadedOutput_[i]->GetPolys()->GetPointer();
    
    float *threadedFiberCoordinates = NULL;
    if(VisibleFibers){
      threadedFiberCoordinates =
        (float *) threadedTextureCoordinates_[i]->GetVoidPointer(0);
    }
    
    int *threadedEdgeIds = NULL;
    if(PolygonSegmentation){
      threadedEdgeIds = (int *) threadedEdgeIds_[i]->GetVoidPointer(0);
    }
      
    for(int j = 0; j < threadedCells_[i]->GetNumberOfCells(); j++){
      
      // number of points in the cell
      cells[4*triangleOffsets[i] + 4*j] = threadedCells[4*j];
     
      for(int k = 0; k < 3; k++){
        cells[4*triangleOffsets[i] + 4*j + 1 + k] = 
          threadedCells[4*j + 1 + k] + 3*triangleOffsets[i];
      }
      
      if(PolygonSegmentation){
        edgeIdsPointer[triangleOffsets[i] + j] = threadedEdgeIds[j];
      }
    }
    
    for(int j = 0; j < threadedPoints_[i]->GetNumberOfPoints(); j++){
      
      for(int k = 0; k < 3; k++){
        points[9*triangleOffsets[i] + 3*j + k] = threadedPoints[3*j + k];
      }
      
      if(VisibleFibers){
        reducedFiberPointer[6*triangleOffsets[i] + 2*j] = 
          threadedFiberCoordinates[2*j];
        reducedFiberPointer[6*triangleOffsets[i] + 2*j + 1] = 
          threadedFiberCoordinates[2*j + 1];
      }
    }
  }
  
  reducedSurface->GetPolys()->SetCells(triangleNumber, cellList);
  if(VisibleFibers){
    reducedSurface->GetPointData()->SetTCoords(reducedFiberCoordinates);
  }
  if(PolygonSegmentation){
    reducedEdgeIds->SetName("Polygon Edge Id");
    reducedSurface->GetCellData()->AddArray(reducedEdgeIds);
  }
  
  if(!Manifold){
    output->ShallowCopy(reducedSurface);
  }
  else{
    vtkSmartPointer<vtkCleanUnstructuredGrid> cleaner = 
      vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfacer = 
      vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
      
    cleaner->SetInputData(reducedSurface);
    cleaner->Update();
    
    surfacer->SetInputData(cleaner->GetOutput());
    surfacer->Update();
    
    output->ShallowCopy(surfacer->GetOutput());
  }
  
  {
    stringstream msg;
    msg << "[vtkFiberSurfaceFilter] Reduce step done in "
      << timer.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), 2);
  }
  
  return 0;
}


int vtkFiberSurfaceFilter::RequestData(vtkInformation *request, 
  vtkInformationVector **in_vec, vtkInformationVector *out_vec){

  if(ThreadNumber <= 0){
    ThreadNumber = 1;
  }
  
  input = vtkUnstructuredGrid::GetData(in_vec[0]);

  vtkUnstructuredGrid *line_cells = vtkUnstructuredGrid::GetData(in_vec[1]);

  output = vtkPolyData::GetData(out_vec);
 
#ifndef withOpenMP
  if(ThreadNumber != 1){
    stringstream msg;
    msg << "[vtkFiberSurfaceFilter]" << endl;
    msg << "[vtkFiberSurfaceFilter] Not compiled with OpenMP support!!!" 
      << endl;
    msg << "[vtkFiberSurfaceFilter] Defaulting back to 1 thread." << endl;
    msg << "[vtkFiberSurfaceFilter]" << endl;
    ThreadNumber = 1;
    dMsg(cerr, msg.str(), 2);
  }
#else  
  {
    stringstream msg;
    msg << "[vtkFiberSurfaceFilter] Using " << ThreadNumber
      << " thread(s)..." << endl;
    dMsg(cout, msg.str(), 2);
  }
#endif

#ifndef withBVH
  if(Implementation == BVH){
    stringstream msg;
    for(int i = 0; i < 10; i++)
      msg << "[vtkFiberSurfaceFilter]" << endl;
    msg << "[vtkFiberSurfaceFilter] Not compiled with BVH support!" << endl;
    msg 
      << "[vtkFiberSurfaceFilter] If your CPU supports AVX features, "
      << "re-run cmake with the option '-DwithBVH=ON'" << endl;
    msg << "[vtkFiberSurfaceFilter] Defaulting back to the octree "
      << "implementation..." << endl;
    for(int i = 0; i < 10; i++)
      msg << "[vtkFiberSurfaceFilter]" << endl;
    Implementation = OCTREE;
    dMsg(cerr, msg.str(), 1);
  }
#endif
 
  // pre-processing
  if(Implementation == OCTREE){
    // octree 
    if(!octree_){
      // retrieve info about the input data-set
      vtkDataArray *fields[2] = { NULL, NULL };

      if(DataUComponent.size()){
        fields[0] = input->GetPointData()->GetArray(DataUComponent.data());
      }
      if(DataVComponent.size()){
        fields[1] = input->GetPointData()->GetArray(DataVComponent.data());
      }
      
      // back up option
      if(!fields[0]){
        fields[0] = input->GetPointData()->GetArray(0);
      }
      if(!fields[1]){
        fields[1] = input->GetPointData()->GetArray(1);
      }
      
      dataType0_ = fields[0]->GetDataType();
      dataType1_ = fields[1]->GetDataType();
      
      switch(dataType0_){
        
        case VTK_CHAR:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildOctree<char, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;

        case VTK_DOUBLE:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildOctree<double, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;
          
        case VTK_FLOAT:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildOctree<float, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;
          
        case VTK_INT:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildOctree<int, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;
          
        case VTK_UNSIGNED_CHAR:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildOctree<unsigned char, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;
          
        default:
          {
            stringstream msg;
            msg << "[vtkFiberSurfaceFilter] Unsupported data type :(" << endl;
            dMsg(cerr, msg.str(), 1);
          }
          break;
      }
    }
  }
  if(Implementation == BVH){
    
    if(!BVH_){
      
      // retrieve info about the input data-set
      vtkDataArray *fields[2] = { NULL, NULL };

      if(DataUComponent.size()){
        fields[0] = input->GetPointData()->GetArray(DataUComponent.data());
      }
      if(DataVComponent.size()){
        fields[1] = input->GetPointData()->GetArray(DataVComponent.data());
      }
      
      // back up option
      if(!fields[0]){
        fields[0] = input->GetPointData()->GetArray(0);
      }
      if(!fields[1]){
        fields[1] = input->GetPointData()->GetArray(1);
      }
      
      dataType0_ = fields[0]->GetDataType();
      dataType1_ = fields[1]->GetDataType();
      
      switch(dataType0_){
        
        case VTK_CHAR:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildBVH<char, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;

        case VTK_DOUBLE:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildBVH<double, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;
          
        case VTK_FLOAT:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildBVH<float, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;
          
        case VTK_INT:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildBVH<int, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;
          
        case VTK_UNSIGNED_CHAR:
          switch(dataType1_){
            vtkTemplateMacro(
              (buildBVH<unsigned char, VTK_TT>(input, 
                fields[0]->GetVoidPointer(0), fields[1]->GetVoidPointer(0)))
            );
          }
          break;
          
        default:
          {
            stringstream msg;
            msg << "[vtkFiberSurfaceFilter] Unsupported data type :(" << endl;
            dMsg(cerr, msg.str(), 1);
          }
          break;
      }
    }
  }
  
  Timer t;
  
  double total_length = 0;
  
  switch(Implementation){
      break;
    case REGULAR: 
      doIt(input, line_cells, output);
      break;
    case OCTREE:
      doItOctree(input, line_cells, output);
      break;
    case BVH:
      doItBVH(input, line_cells, output);
      break;
  }
  
  {
    stringstream msg;
    msg << "[vtkFiberSurfaceFilter] Fiber surface computed in "
      << t.getElapsedTime() << " s. TOTAL ("
      << output->GetNumberOfPoints() << " #v, "
      << output->GetNumberOfCells() << " #f)" << endl;
    dMsg(cout, msg.str(), 2);
  }
  
  // post-process from now on
  if(VisibleFibers){
    vtkDataArray *t_coords_array = output->GetPointData()->GetTCoords();
    // julien: mandatory, otherwise ParaView doesn't see it.
    t_coords_array->SetName("Fiber Coordinates");
    const vtkIdType t_coords_size = t_coords_array->GetNumberOfTuples();
    for (vtkIdType i = 0; i != t_coords_size; ++i)
      t_coords_array->SetComponent(i, 1, 
        t_coords_array->GetComponent(i, 1)/total_length);
  }

  return 1;
}


// PRIVATE
const char vtkFiberSurfaceFilter::name[] = "vtkFiberSurfaceFilter";
