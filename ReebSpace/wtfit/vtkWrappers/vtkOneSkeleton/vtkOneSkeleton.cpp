#include                  <vtkOneSkeleton.h>

vtkOneSkeleton::vtkOneSkeleton(){

  debugLevel_ = 0;
}

int vtkOneSkeleton::buildEdgeList(vtkDataSet *input, 
  vector<pair<int, int> > &edgeList,
  const bool &isTriangulation) const{
    
#ifndef withKamikaze
  if(!input)
    return -1;
#endif
    
  if(isTriangulation){
    if(input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
      return 
        buildTriangulationEdgeList((vtkUnstructuredGrid *) input, edgeList);
    }
    if(input->GetDataObjectType() == VTK_POLY_DATA){
      return buildTriangulationEdgeList((vtkPolyData *) input, edgeList);
    }
  }
  
  Timer t;
  
  vector<vector<vector<int> > > threadedEdgeTable(threadNumber_);
  vector<vtkSmartPointer<vtkIdList> > threadedStarList(threadNumber_);
  vector<vtkSmartPointer<vtkGenericCell> > threadedCells(threadNumber_);
  
  for(int i = 0; i < (int) threadedEdgeTable.size(); i++){
    threadedEdgeTable[i].resize(input->GetNumberOfPoints());
    threadedStarList[i] = vtkSmartPointer<vtkIdList>::New();
    threadedStarList[i]->Allocate(8);
    threadedCells[i] = vtkSmartPointer<vtkGenericCell>::New();
  }
  
  // make these functions thread-safe
  input->GetPointCells(0, threadedStarList[0]);
  input->GetCell(0, threadedCells[0]);
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < input->GetNumberOfPoints(); i++){
    
    int threadId = 0;
    
#ifdef withOpenMP
    threadId = omp_get_thread_num();
#endif
    
    input->GetPointCells(i, threadedStarList[threadId]);
    
    for(int j = 0; j < threadedStarList[threadId]->GetNumberOfIds(); j++){
      
      input->GetCell(threadedStarList[threadId]->GetId(j), 
        threadedCells[threadId]);
      
      if(threadedCells[threadId]->GetNumberOfEdges()){
        
        for(int k = 0; k < threadedCells[threadId]->GetNumberOfEdges(); k++){
          
          vtkCell *edge = threadedCells[threadId]->GetEdge(k);
          
          pair<int, int> edgeIds;
          edgeIds.first = edge->GetPointId(0);
          edgeIds.second = edge->GetPointId(1);
          
          // canonical layout
          if(edgeIds.second < edgeIds.first){
            edgeIds.first = edgeIds.second;
            edgeIds.second = edge->GetPointId(0);
          }
          
             
          
          
          // search if it already exists
          bool hasFound = false;
          for(int l = 0; 
            l < (int) threadedEdgeTable[threadId][edgeIds.first].size(); l++){
            
            if(edgeIds.second == threadedEdgeTable[threadId][edgeIds.first][l]){
              
              hasFound = true;
              break;
            }
          }
          if(!hasFound){
            threadedEdgeTable[
              threadId][edgeIds.first].push_back(edgeIds.second);
          }
        }
      }
      else{
        // the input is a 1-mesh (no edge within an edge)
        pair<int, int> edgeIds;
        edgeIds.first = threadedCells[threadId]->GetPointId(0);
        edgeIds.second = threadedCells[threadId]->GetPointId(1);
        
        // canonical layout
        if(edgeIds.second < edgeIds.first){
          edgeIds.first = edgeIds.second;
          edgeIds.second = threadedCells[threadId]->GetPointId(0);
        }
        
        // search if it already exists
        bool hasFound = false;
        for(int l = 0; 
          l < (int) threadedEdgeTable[threadId][edgeIds.first].size(); l++){
          
          if(edgeIds.second == threadedEdgeTable[threadId][edgeIds.first][l]){
            hasFound = true;
            break;
          }
        }
        if(!hasFound){
          threadedEdgeTable[
            threadId][edgeIds.first].push_back(edgeIds.second);
        }
      }
    }
  }
  
  // now merge the thing
  int edgeCount = 0;
  vector<vector<int> > edgeTable(input->GetNumberOfPoints());
  for(int i = 0; i < (int) threadedEdgeTable.size(); i++){
    
    for(int j = 0; j < (int) threadedEdgeTable[i].size(); j++){
      
      for(int k = 0; k < (int) threadedEdgeTable[i][j].size(); k++){
        
        // search if it already exists
        bool hasFound = false;
        
        for(int l = 0; l < (int) edgeTable[j].size(); l++){
          if(edgeTable[j][l] == threadedEdgeTable[i][j][k]){
            hasFound = true;
            break;
          }
        }
        if(!hasFound){
          edgeTable[j].push_back(threadedEdgeTable[i][j][k]);
          edgeCount++;
        }
      }
    }
  }
 
  edgeList.resize(edgeCount);
  edgeCount = 0;
  for(int i = 0; i < (int) edgeTable.size(); i++){
    
    for(int j = 0; j < (int) edgeTable[i].size(); j++){
      
      edgeList[edgeCount].first = i;
      edgeList[edgeCount].second = edgeTable[i][j];
      edgeCount++;
    }
  }
  
  {
    stringstream msg;
    msg << "[vtkOneSkeleton] Edge-list built in " 
      << t.getElapsedTime() << " s. (" << edgeList.size()
      << " edges, (" << threadNumber_
      << " thread(s))" 
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 159.781 s
  // 24 threads:  34.2091 s [not efficient in parallel]
  
  return 0;
}

int vtkOneSkeleton::buildEdgeStars(vtkDataSet* input, 
  vector<vector<int> > &starList,
  vector<pair<int, int> > *edgeList,
  vector<vector<int> > *vertexStars,
  const bool &isTriangulation) const{

#ifndef withKamikaze
  if(!input)
    return -1;
#endif
    
  if(isTriangulation){
    
    if(input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
      return 
        buildTriangulationEdgeStars((vtkUnstructuredGrid *) input, starList,
          edgeList, vertexStars);
    }
    if(input->GetDataObjectType() == VTK_POLY_DATA){
      return buildTriangulationEdgeStars((vtkPolyData *) input, starList,
        edgeList, vertexStars);
    }
  }
  
  Timer t;
  
  bool localEdgeListAlloc = false;
  vector<pair<int, int> > *localEdgeList = edgeList;
  
  if(!localEdgeList){
    localEdgeList = new vector<pair<int, int> >();
    localEdgeListAlloc = true;
  }
  
  if(!localEdgeList->size()){
    buildEdgeList(input, *localEdgeList);
  }
  
  starList.resize(localEdgeList->size());
 
  vector<vtkSmartPointer<vtkIdList> > vertex0starList, vertex1starList;
  
  vertex0starList.resize(threadNumber_);
  vertex1starList.resize(threadNumber_);
  
  for(int i = 0; i < (int) vertex0starList.size(); i++){
    vertex0starList[i] = vtkSmartPointer<vtkIdList>::New();
    vertex1starList[i] = vtkSmartPointer<vtkIdList>::New();
  }
  
  // make the call threadsafe
  input->GetPointCells(0, vertex0starList[0]);
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < (int) localEdgeList->size(); i++){
    
    int threadId = 0;
    
#ifdef withOpenMP
    threadId = omp_get_thread_num();
#endif
    
    input->GetPointCells(
      (*localEdgeList)[i].first, vertex0starList[threadId]);
    input->GetPointCells(
      (*localEdgeList)[i].second, vertex1starList[threadId]);
    
    // merge the two
    for(int j = 0; j < vertex0starList[threadId]->GetNumberOfIds(); j++){
      
      bool hasFound = false;
      for(int k = 0; k < vertex1starList[threadId]->GetNumberOfIds(); k++){
        if(vertex0starList[threadId]->GetId(j) 
          == vertex1starList[threadId]->GetId(k)){
        
          hasFound = true;
          break;
        }
      }
      if(hasFound){
        starList[i].push_back(vertex0starList[threadId]->GetId(j));
      }
    }
  }
  
  if(localEdgeListAlloc)
    delete localEdgeList;
  
  {
    stringstream msg;
    msg << "[vtkOneSkeleton] List of edge stars built in " 
      << t.getElapsedTime() << " s. (" << starList.size()
      << " edges, (" << threadNumber_
      << " thread(s))"<< endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 72.1875 s
  // 24 threads: 23.0986 s (alright) 
  
  return 0;
}
