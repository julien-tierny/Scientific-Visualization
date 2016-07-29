#include                  <vtkZeroSkeleton.h>

vtkZeroSkeleton::vtkZeroSkeleton(){

  debugLevel_ = 0;
}

int vtkZeroSkeleton::buildVertexEdges(vtkDataSet *input, 
  vector<vector<int> > &vertexEdges, vector<pair<int,int> > *edgeList, 
  const bool &isTriangulation) const{

#ifndef withKamikaze
  if(!input)
    return -1;
#endif
  
  bool localEdgeListAlloc = false;
  vector<pair<int, int> > *localEdgeList = edgeList;
  
  if(!localEdgeList){
    localEdgeList = new vector<pair<int, int> >();
    localEdgeListAlloc = true;
  }
  
  if(!localEdgeList->size()){
    vtkOneSkeleton oneSkeleton;
    oneSkeleton.setWrapper(this);
    oneSkeleton.buildEdgeList(input, *localEdgeList, isTriangulation);
  }
  
  ZeroSkeleton zeroSkeleton;
  zeroSkeleton.setWrapper(this);
  zeroSkeleton.buildVertexEdges(input->GetNumberOfPoints(), 
    *localEdgeList, vertexEdges);
  
  if(localEdgeListAlloc)
    delete localEdgeList;
  
  return 0;
}


int vtkZeroSkeleton::buildVertexLinks(vtkDataSet *input, 
  vector<vector<long long int > > &vertexLinks,
  vector<vector<int> > *vertexStars,
  const bool &isTriangulation) const{

#ifndef withKamikaze
  if(!input)
    return -1;
#endif
  
  if(isTriangulation){
    if(input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
      return 
        buildTriangulationVertexLinks((vtkUnstructuredGrid *) input, 
          vertexLinks);
    }
    if(input->GetDataObjectType() == VTK_POLY_DATA){
      return buildTriangulationVertexLinks((vtkPolyData *) input, 
        vertexLinks);
    }
  }
  
  Timer t;
 
  // init
  vertexLinks.resize(input->GetNumberOfPoints());
 
  vector<vtkIdList *> starIdsVector(threadNumber_);
  vector<vtkGenericCell *> cellVector(threadNumber_);
  vector<vector<long long int> > faceIds(threadNumber_);
 
  for(int i = 0; i < threadNumber_; i++){
    starIdsVector[i] = vtkIdList::New();
    cellVector[i] = vtkGenericCell::New();
    faceIds[i].resize(4);
    faceIds[i][0] = 3;
  }
  
  // make these functions thread-safe
  input->GetPointCells(0, starIdsVector[0]);
  input->GetCell(0, cellVector[0]);
 
  // init 
  for(int i = 0; i < input->GetNumberOfPoints(); i++){
    input->GetPointCells(i, starIdsVector[0]);
    vertexLinks[i].reserve(
      starIdsVector[0]->GetNumberOfIds()*4);
  }
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < input->GetNumberOfPoints(); i++){
    
    int threadId = 0;
    
#ifdef withOpenMP
    threadId = omp_get_thread_num();
#endif
    
    // this call is implemented differently for regular grids and 
    // unstructured grids. no way to get a pointer easily.
    input->GetPointCells(i, starIdsVector[threadId]);
    
    // in 3D, let's loop over the tets
    // pre-allocate the memory
    for(int j = 0; j < starIdsVector[threadId]->GetNumberOfIds(); j++){
      input->GetCell(starIdsVector[threadId]->GetId(j), cellVector[threadId]);
      
      // in 3D, let's loop over the triangles of the tets
      switch(cellVector[threadId]->GetCellDimension()){
        
        case 3:
          {
            for(int k = 0; k < cellVector[threadId]->GetNumberOfFaces(); k++){
              vtkCell *face = cellVector[threadId]->GetFace(k);
              
              bool hasVertex = false;
              for(int l = 0; l < face->GetNumberOfPoints(); l++){
                int neighborId = face->GetPointId(l);
                if(neighborId == i){
                  hasVertex = true;
                  break;
                }
                faceIds[threadId][l + 1] = neighborId;
              }
              if(!hasVertex){
                for(int l = 0; l < (int) faceIds[threadId].size(); l++){
                  vertexLinks[i].push_back(faceIds[threadId][l]);
                }
//                 vertexLinks[i][j] = faceIds[threadId];
                break;
              }
            }
          } 
          break;
          
        case 2:
          {
            for(int k = 0; k < cellVector[threadId]->GetNumberOfEdges(); k++){
              vtkCell *edge = cellVector[threadId]->GetEdge(k);
              
              bool hasVertex = false;
              if(faceIds[threadId].size() == 4){
                faceIds[threadId].resize(3);
                faceIds[threadId][0] = 2;
              }
              for(int l = 0; l < edge->GetNumberOfPoints(); l++){
                int neighborId = edge->GetPointId(l);
                if(neighborId == i){
                  hasVertex = true;
                  break;
                }
                faceIds[threadId][l + 1] = neighborId;
              }
              if(!hasVertex){
                for(int l = 0; l < (int) faceIds[threadId].size(); l++){
                  vertexLinks[i].push_back(faceIds[threadId][l]);
                }
//                 vertexLinks[i][j] = faceIds[threadId];
                break;
              }
            }
          }
          break;
          
        case 1:
          {
            for(int k = 0; k < cellVector[threadId]->GetNumberOfPoints(); k++){
              int neighborId = cellVector[threadId]->GetPointId(k);
              if(faceIds[threadId].size() == 3){
                faceIds[threadId].resize(2);
                faceIds[threadId][0] = 1;
              }
              if(neighborId != i){
                faceIds[threadId][1] = neighborId;
                for(int l = 0; l < (int) faceIds[threadId].size(); l++){
                  vertexLinks[i].push_back(faceIds[threadId][l]);
                }
//                 vertexLinks[i][j] = faceIds[threadId];
                break;
              }
            }
          }
          break;
      }
    }
  }
  
  for(int i = 0; i < threadNumber_; i++){
    cellVector[i]->Delete();
    starIdsVector[i]->Delete();
  }
  
  {
    stringstream msg;
    msg << "[vtkZeroSkeleton] Vertex links built in " 
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 87.28 s
  // 24 threads: 22.65 s
  
  return 0;
}

int vtkZeroSkeleton::buildVertexNeighbors(vtkDataSet *input, 
  vector<vector<int> > &vertexNeighbors,
  const bool &isTriangulation) const{ 

#ifndef withKamikaze
  if(!input)
    return -1;
#endif
    
  if(isTriangulation){
    if(input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
      return 
        buildTriangulationVertexNeighbors((vtkUnstructuredGrid *) input, 
          vertexNeighbors);
    }
    if(input->GetDataObjectType() == VTK_POLY_DATA){
      return buildTriangulationVertexNeighbors((vtkPolyData *) input, 
        vertexNeighbors);
    }
  }
  
  Timer t;
 
  // getting one skeleton
  vertexNeighbors.resize(input->GetNumberOfPoints());
  vector<vtkIdList * > starIdsVector(threadNumber_);
  vector<vtkGenericCell *> cellVector(threadNumber_);
  
  
  for(int i = 0; i < (int) starIdsVector.size(); i++){
    starIdsVector[i] = vtkIdList::New();
    starIdsVector[i]->Allocate(32);
  }
  for(int i = 0; i < (int) cellVector.size(); i++){
    cellVector[i] = vtkGenericCell::New();
  }
  
  // pre-allocate
  for(int i = 0; i < (int) vertexNeighbors.size(); i++){
    vertexNeighbors[i].reserve(32);
  }
  
  // make these functions thread-safe
  input->GetPointCells(0, starIdsVector[0]);
  input->GetCell(0, cellVector[0]);

#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < (int) vertexNeighbors.size(); i++){
    
    int threadId = 0;
    
#ifdef withOpenMP
    threadId = omp_get_thread_num();
#endif
    
    input->GetPointCells(i, starIdsVector[threadId]);
    
    for(int j = 0; j < starIdsVector[threadId]->GetNumberOfIds(); j++){
      
      // this call is implemented differently for regular grids and 
      // unstructured grids. no way to get a pointer easily.
      input->GetCell(starIdsVector[threadId]->GetId(j), cellVector[threadId]);
      
      if(cellVector[threadId]->GetNumberOfEdges()){
        for(int k = 0; k < cellVector[threadId]->GetNumberOfEdges(); k++){
          vtkCell *edge = cellVector[threadId]->GetEdge(k);
          
          int vertexId0 = edge->GetPointId(0);
          int vertexId1 = edge->GetPointId(1);
          
          if((vertexId0 == i)||(vertexId1 == i)){
            if(vertexId0 == i){
              vertexNeighbors[i].push_back(vertexId1);
            }
            else{
              vertexNeighbors[i].push_back(vertexId0);
            }
          }
        }
      }
      else{
        // the input is a 1-mesh, (so there's no edge within an edge)
        for(int k = 0; k < cellVector[threadId]->GetNumberOfPoints(); k++){
          int vertexId0 = cellVector[threadId]->GetPointId(0);
          int vertexId1 = cellVector[threadId]->GetPointId(1);
          
          if(vertexId0 == i){
            vertexNeighbors[i].push_back(vertexId1);
          }
          else{
            vertexNeighbors[i].push_back(vertexId0);
          }
        }
      }
    }
    
    // remove duplicates
    vector<int>::iterator it;
    sort(vertexNeighbors[i].begin(), vertexNeighbors[i].end());
    it = unique(vertexNeighbors[i].begin(), vertexNeighbors[i].end());
    vertexNeighbors[i].resize(distance(vertexNeighbors[i].begin(), it));
  }
 
  for(int i = 0; i < (int) starIdsVector.size(); i++){
    starIdsVector[i]->Delete();
  }
  for(int i = 0; i < (int) cellVector.size(); i++){
    cellVector[i]->Delete();
  }
  
  {
    stringstream msg;
    msg << "[vtkZeroSkeleton] Vertex neighbors built in " 
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
 
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 112.165 s
  // 24 threads: 18.66 s (~ x6)
 
  return 0;
}

int vtkZeroSkeleton::buildVertexStars(vtkDataSet *input, 
  vector<vector<int > > &vertexStars, const bool &isTriangulation) const{ 

#ifndef withKamikaze
  if(!input)
    return -1;
#endif
    
  if(isTriangulation){
    if(input->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
      return 
        buildTriangulationVertexStars((vtkUnstructuredGrid *) input, 
          vertexStars);
    }
    if(input->GetDataObjectType() == VTK_POLY_DATA){
      return buildTriangulationVertexStars((vtkPolyData *) input, 
        vertexStars);
    }
  }
  
  Timer t;
 
  // init
  vertexStars.resize(input->GetNumberOfPoints());
 
  vector<vtkIdList *> starIdsVector(threadNumber_);
 
  for(int i = 0; i < threadNumber_; i++){
    starIdsVector[i] = vtkIdList::New();
    starIdsVector[i]->Allocate(16);
  }

  // make the call thread safe
  input->GetPointCells(0, starIdsVector[0]);
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < input->GetNumberOfPoints(); i++){
    
    int threadId = 0;
#ifdef withOpenMP
    threadId = omp_get_thread_num();
#endif
    
    input->GetPointCells(i, starIdsVector[threadId]);
    vertexStars[i].resize(starIdsVector[threadId]->GetNumberOfIds());
    
    for(int j = 0; j < (int) vertexStars[i].size(); j++)
      vertexStars[i][j] = starIdsVector[threadId]->GetId(j);
  }
  
  for(int i = 0; i < threadNumber_; i++){
    starIdsVector[i]->Delete();
  }
  
  {
    stringstream msg;
    msg << "[vtkZeroSkeleton] Vertex stars built in " 
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  // ethaneDiol.vtu, 8.7Mtets, hal9000 (12coresHT)
  // 1 thread: 1.44 s
  // 24 threads: 1.57 s 
  
  // ethaneDiolMedium.vtu, 70Mtets, hal9000 (12coresHT)
  // 1 thread: 14.68 s
  // 24 threads: 12.34 s 
  
  return 0;
}