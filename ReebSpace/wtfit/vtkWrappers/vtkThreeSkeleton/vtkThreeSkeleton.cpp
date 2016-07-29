#include                  <vtkThreeSkeleton.h>

vtkStandardNewMacro(vtkThreeSkeleton)

vtkThreeSkeleton::vtkThreeSkeleton(){

  // init
  debugLevel_ = 0;
}

vtkThreeSkeleton::~vtkThreeSkeleton(){

}

int vtkThreeSkeleton::buildCellEdges(vtkDataSet *input, 
  vector<vector<int> > &cellEdges, vector<pair<int, int> > *edgeList, 
  vector<vector<int> > *vertexEdges, const bool &isTriangulation) const {

#ifndef withKamikaze
  if(!input)
    return -1;
  if(!input->GetNumberOfCells())
    return -2;
#endif
    
  Timer t;
  
  bool localEdgeListAlloc = false, localVertexEdgesAlloc = false;
  vector<pair<int, int> > *localEdgeList = edgeList;
  vector<vector<int> > *localVertexEdges = vertexEdges;
  
  if(!localEdgeList){
    localEdgeList = new vector<pair<int, int> >();
    localEdgeListAlloc = true;
  }
  
  if(!localEdgeList->size()){
    
    vtkOneSkeleton oneSkeleton;
    oneSkeleton.setWrapper(this);
    oneSkeleton.buildEdgeList(input, *localEdgeList, isTriangulation);
  }
  
  if(!localVertexEdges){
    localVertexEdges = new vector<vector<int> >();
    localVertexEdgesAlloc = true;
  }
  
  if(!localVertexEdges->size()){
    
    vtkZeroSkeleton zeroSkeleton;
    zeroSkeleton.setWrapper(this);
    zeroSkeleton.buildVertexEdges(input, 
      *localVertexEdges, localEdgeList, isTriangulation);
  }
  
  // TODO: fast implementation for triangulations (using edge stars)
  
  vector<vtkSmartPointer<vtkGenericCell> > threadedCells(threadNumber_);
  for(int i = 0; i < threadNumber_; i++){
    threadedCells[i] = vtkSmartPointer<vtkGenericCell>::New();
  }
  
  cellEdges.resize(input->GetNumberOfCells());
  for(int i = 0; i < (int) cellEdges.size(); i++){
    // pre-reserve based on tet-meshes
    cellEdges[i].reserve(6);
  }
  
  {
    // make vtk calls thread safe
    input->GetCell(0, threadedCells[0]);
  }

 
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < input->GetNumberOfCells(); i++){
    
    int threadId = 0;
#ifdef withOpenMP
    threadId = omp_get_thread_num();
#endif
    
    input->GetCell(i, threadedCells[threadId]);
    
    for(int j = 0; j < threadedCells[threadId]->GetNumberOfPoints(); j++){
      
      for(int k = j + 1; k < threadedCells[threadId]->GetNumberOfPoints(); k++){
      
        int vertexId0 = threadedCells[threadId]->GetPointId(j);
        int vertexId1 = threadedCells[threadId]->GetPointId(k);
        
        int edgeId = -1;
        for(int l = 0; l < (int) (*localVertexEdges)[vertexId0].size(); l++){
          
          for(int m = 0; 
            m < (int) (*localVertexEdges)[vertexId1].size(); m++){
            
            if((*localVertexEdges)[vertexId0][l] == 
              (*localVertexEdges)[vertexId1][m]){
              edgeId = (*localVertexEdges)[vertexId0][l];
              break;
            }
          }
          if(edgeId != -1)
            break;
        }
        cellEdges[i].push_back(edgeId);
      }
    }
  }
  
  if(localEdgeListAlloc)
    delete localEdgeList;
  if(localVertexEdgesAlloc)
    delete localVertexEdges;
    
  {
    stringstream msg;
    msg << "[vtkThreeSkeleton] Cell edges built in " 
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // (including memory allocation, ~ 5s.)
  // 1 thread: 92.0522 s
  // 24 threads: 11.1999 s
  
  return 0;
}


int vtkThreeSkeleton::buildCellNeighbors(vtkDataSet *input, 
  vector<vector<int> > &cellNeighbors, 
  vector<vector<int> > *vertexStars,
  vector<vector<int> > *triangleStars,
  const bool &isTriangulation) const {

#ifndef withKamikaze
  if(!input)
    return -1;
  if(!input->GetNumberOfCells())
    return -2;
#endif
  
  if(isTriangulation){
    if(input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
      return 
        buildTriangulationCellNeighbors(
          (vtkUnstructuredGrid *) input, 
            cellNeighbors, vertexStars, triangleStars);
    }
    if(input->GetDataObjectType() == VTK_POLY_DATA){
      return 
        buildTriangulationCellNeighbors(
          (vtkPolyData *) input, 
            cellNeighbors, vertexStars, triangleStars);
    }
  }
  
  Timer t;
  
  int domainDimension = 
    (input->GetCell(0)->GetNumberOfFaces() > 0 ? 3
    : (input->GetCell(0)->GetNumberOfEdges() > 0 ? 2
      : 1));
      
  cellNeighbors.resize(input->GetNumberOfCells());
  
  if(domainDimension == 3){
    for(int i = 0; i < (int) cellNeighbors.size(); i++)
      cellNeighbors[i].reserve(input->GetCell(i)->GetNumberOfFaces());
  }
  else if(domainDimension == 2){
    for(int i = 0; i < (int) cellNeighbors.size(); i++)
      cellNeighbors[i].reserve(input->GetCell(i)->GetNumberOfEdges());
  }
  else if(domainDimension == 1){
    for(int i = 0; i < (int) cellNeighbors.size(); i++)
      cellNeighbors[i].reserve(input->GetCell(i)->GetNumberOfPoints());
  }
  // NOTE: reserve is the best we can do (boundary cells...)
 
  vector<vtkSmartPointer<vtkIdList> > threadedIdList(threadNumber_);
  for(int i = 0; i < (int) threadedIdList.size(); i++)
    threadedIdList[i] = vtkSmartPointer<vtkIdList>::New();
  
  vector<vtkSmartPointer<vtkGenericCell> > threadedCells(threadNumber_);
  for(int i = 0; i < (int) threadedCells.size(); i++)
    threadedCells[i] = vtkSmartPointer<vtkGenericCell>::New();
  
  // make these calls thread-safe
  {
    input->GetCell(0, threadedCells[0]);
    vtkCell *face = threadedCells[0]->GetFace(0);
    
    if(face){
      vtkIdList *vertexIds = face->GetPointIds();
      input->GetCellNeighbors(0, vertexIds, threadedIdList[0]);
    }
    else{
      vtkCell *edge = threadedCells[0]->GetEdge(0);
      vtkIdList *vertexIds = edge->GetPointIds();
      input->GetCellNeighbors(0, vertexIds, threadedIdList[0]);
    }
  }
  
  if(domainDimension == 3){
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < (int) cellNeighbors.size(); i++){
    
      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
#endif
      
      input->GetCell(i, threadedCells[threadId]);
      vtkIdList *vertexIds = NULL;
    
      for(int j = 0; j < threadedCells[threadId]->GetNumberOfFaces(); j++){
        
        vtkCell *face = threadedCells[threadId]->GetFace(j);
        vertexIds = face->GetPointIds();
        
        input->GetCellNeighbors(i, vertexIds, threadedIdList[threadId]);
        
        if(threadedIdList[threadId]->GetNumberOfIds()){
          cellNeighbors[i].push_back(threadedIdList[threadId]->GetId(0));
        }
      }
    }
  }
  else if(domainDimension == 2){
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < (int) cellNeighbors.size(); i++){
    
      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
#endif
      
      input->GetCell(i, threadedCells[threadId]);
      vtkIdList *vertexIds = NULL;
    
      for(int j = 0; j < threadedCells[threadId]->GetNumberOfEdges(); j++){
        
        vtkCell *edge = threadedCells[threadId]->GetEdge(j);
        vertexIds = edge->GetPointIds();
        
        input->GetCellNeighbors(i, vertexIds, threadedIdList[threadId]);
        
        if(threadedIdList[threadId]->GetNumberOfIds()){
          cellNeighbors[i].push_back(threadedIdList[threadId]->GetId(0));
        }
      }
    }
  }
  else{
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(int i = 0; i < (int) cellNeighbors.size(); i++){
    
      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
#endif
      
      input->GetCell(i, threadedCells[threadId]);
      vtkIdList *vertexIds = vtkIdList::New();
      vertexIds->SetNumberOfIds(1);
    
      for(int j = 0; j < threadedCells[threadId]->GetNumberOfPoints(); j++){
        
        vertexIds->SetId(0, threadedCells[threadId]->GetPointId(j));
        
        input->GetCellNeighbors(i, vertexIds, threadedIdList[threadId]);
        
        if(threadedIdList[threadId]->GetNumberOfIds()){
          cellNeighbors[i].push_back(threadedIdList[threadId]->GetId(0));
        }
      }
      
      vertexIds->Delete();
    }   
  }
  
  {
    stringstream msg;
    msg << "[vtkThreeSkeleton] Cell neighbors built in " 
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 96.2044 s
  // 8 threads: 35.7842 s
  // 24 threads: 27.2799 s (speedup ~3.5)
    
  // ethaneDiol.vtu, 8.7Mtets, hal9000 (12coresHT)
  // 1 thread: 11.3164 s
  // 8 threads: 4.04992 s
  // 24 threads: 3.1558 s (speedup ~3.5)
    
  // ethaneDiol.vtu, 8.7Mtets, vger (4coresHT)
  // 1 thread: 15.7514 s
  // 4 threads: 8.16998 s
    
  return 0;
}


// transmit abort signals -- to copy paste in other wrappers
bool vtkThreeSkeleton::needsToAbort(){
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int vtkThreeSkeleton::updateProgress(const float &progress){

  {
    stringstream msg;
    msg << "[vtkThreeSkeleton] " << progress*100 
      << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  UpdateProgress(progress);
  return 0;
}

int vtkThreeSkeleton::doIt(vtkDataSet *input, vtkDataSet *output){

  return 0;
}

// to adapt if your wrapper does not inherit from vtkUnstructuredGridAlgorithm
int vtkThreeSkeleton::RequestData(vtkInformation *request, 
  vtkInformationVector **inputVector, vtkInformationVector *outputVector){

  Memory m;
  
  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);
  
  doIt(input, output);
  
  {
    stringstream msg;
    msg << "[vtkThreeSkeleton] Memory usage: " << m.getElapsedUsage() 
      << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }
  
  return 1;
}