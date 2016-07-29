#include                  <JacobiSet.h>

template <class dataTypeU, class dataTypeV> 
  JacobiSet<dataTypeU, dataTypeV>::JacobiSet(){

  vertexNumber_ = 0;
  
  uField_ = NULL;
  vField_ = NULL;
  
  tetList_ = NULL;
  
  edgeList_ = NULL;
  edgeFanLinkEdgeLists_ = NULL;
  edgeFans_ = NULL;
  sosOffsets_ = NULL;
  
}

template <class dataTypeU, class dataTypeV> 
  JacobiSet<dataTypeU, dataTypeV>::~JacobiSet(){
  
}

template <class dataTypeU, class dataTypeV> 
  int JacobiSet<dataTypeU, dataTypeV>::connectivityPreprocessing(
    const vector<vector<int> > &edgeStarList,
    vector<vector<pair<int, int> > > &edgeFanLinkEdgeLists,
    vector<vector<long long int> > &edgeFans,
    vector<int> &sosOffsets) const{

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef withKamikaze
  if(!vertexNumber_)
    return -1;
  if(!tetList_)
    return -2;
  if(!edgeList_)
    return -3;
  if(edgeStarList.size() != edgeList_->size())
    return -4;
#endif

  edgeFanLinkEdgeLists.resize(edgeList_->size());
  
  int count = 0;

  // edge triangle fans [each thread writes in a different spot]
  //    for each edge
  //      for each triangle 
  //        list of vertices
  edgeFans.resize(edgeList_->size());
  // pre-allocate memory
  for(int i = 0; i < (int) edgeFans.size(); i++){
    // we store 4 integers per triangle per edge
    edgeFans[i].resize(edgeStarList[i].size()*4);
  }
 
  if(!sosOffsets.size()){
    sosOffsets.resize(vertexNumber_);
    for(int i = 0; i < vertexNumber_; i++){
      sosOffsets[i] = i;
    }
  }
  
  vector<ZeroSkeleton> threadedLinkers(threadNumber_);
  vector<vector<long long int> > threadedLinks(threadNumber_);
  for(int i = 0; i < threadNumber_; i++){
    threadedLinkers[i].setDebugLevel(debugLevel_);
    threadedLinkers[i].setThreadNumber(1);
  }
  
  vector<OneSkeleton> threadedEdgeListers(threadNumber_);
  for(int i = 0; i < threadNumber_; i++){
    threadedEdgeListers[i].setDebugLevel(debugLevel_);
    threadedEdgeListers[i].setThreadNumber(1);
  }
 
  vector<vector<int> > threadedTriangleIds(threadNumber_);
  for(int i = 0; i < (int) threadedTriangleIds.size(); i++){
    threadedTriangleIds[i].resize(4);
    threadedTriangleIds[i][0] = 3;
  }
  
#ifdef withOpenMP
  omp_lock_t writeLock;
  omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(int i = 0; i < (int) edgeList_->size(); i++){

    // avoid any processing if the abort signal is sent
    if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){

      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
#endif
      
      // processing here!
      int pivotVertexId = (*edgeList_)[i].first;
      int otherExtremityId = (*edgeList_)[i].second;
      
      // A) compute triangle fans
      // format: #vertices, id0, id1, id2, etc.
      for(int j = 0; j < (int) edgeStarList[i].size(); j++){
        
        int tetId = edgeStarList[i][j];
        
        // loop over the tet's triangles and add that to the list
        // no need for check, only one triangle verifies this and two tets
        // can't add the same triangle
        for(int k = 0; k < 4; k++){
          
          bool hasPivotVertex = false;
          bool hasOtherExtremity = false;
          for(int l = 0; l < 3; l++){
            threadedTriangleIds[threadId][l + 1] 
              = tetList_[5*tetId + 1 + (l + k)%4];
            if(threadedTriangleIds[threadId][l + 1] == pivotVertexId){
              hasPivotVertex = true;
            }
            if(threadedTriangleIds[threadId][l + 1] == otherExtremityId){
              hasOtherExtremity = true;
            }
          }
          
          if((hasPivotVertex)&&(!hasOtherExtremity)){
            for(int l = 0; l < (int) threadedTriangleIds[threadId].size(); l++){
              edgeFans[i][j*4 + l] = 
                threadedTriangleIds[threadId][l];
            }
            break;
          }
        }
      }
      
      // set-up the link of the edge fan
      threadedLinkers[threadId].buildVertexLink(
        pivotVertexId, edgeFans[i].size()/4, 
          edgeFans[i].data(), threadedLinks[threadId]);
      
      // now compute the edge list of the link
      threadedEdgeListers[threadId].buildEdgeSubList(
        edgeFans[i].size()/4, edgeFans[i].data(),
        edgeFanLinkEdgeLists[i]);
      
      // update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > advancedInfoMsg){
#ifdef withOpenMP
        omp_set_lock(&writeLock);
#endif
        if((wrapper_)
          &&(!(count % ((vertexNumber_)/10)))){
          wrapper_->updateProgress((count + 1.0)
            /vertexNumber_);
        }

        count++;
#ifdef withOpenMP
        omp_unset_lock(&writeLock);
#endif
      }
    }
  }
   
#ifdef withOpenMP
  omp_destroy_lock(&writeLock);
#endif
  
  {
    stringstream msg;
    msg << "[JacobiSet] Edge-fans computed in "
      << t.getElapsedTime() << " s. ("
      << edgeList_->size()
      << " edges)" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }
  
  return 0;
}

template <class dataTypeU, class dataTypeV> 
  int JacobiSet<dataTypeU, dataTypeV>::execute(
    vector<pair<int, char> > &jacobiSet){

  Timer t;
  
  // check the consistency of the variables -- to adapt
#ifndef withKamikaze
  if(!vertexNumber_)
    return -1;
  if(!uField_)
    return -2;
  if(!vField_)
    return -3;
  if(!edgeList_)
    return -4;
  if(!edgeFanLinkEdgeLists_)
    return -5;
  if(!edgeFans_)
    return -6;
  if(!sosOffsets_)
    return -7;
#endif

  int count = 0;
  
  jacobiSet.clear();
  
  dataTypeU *uField = (dataTypeU *) uField_;
  dataTypeV *vField = (dataTypeV *) vField_;
 
  // distance fields (not really memory efficient)
  // for each thread
  //      for each vertex: distance field map
  vector<vector<double> > threadedDistanceField(threadNumber_);
  for(int i = 0; i < (int) threadedDistanceField.size(); i++){
    threadedDistanceField[i].resize(vertexNumber_);
  }
  
  vector<ScalarFieldCriticalPoints<double> > 
    threadedCriticalPoints(threadNumber_);
  for(int i = 0; i < threadNumber_; i++){
    threadedCriticalPoints[i].setDomainDimension(2);
    threadedCriticalPoints[i].setScalarValues(
      threadedDistanceField[i].data());
    threadedCriticalPoints[i].setVertexNumber(vertexNumber_);
    threadedCriticalPoints[i].setSosOffsets(sosOffsets_);
  }

  vector<vector<pair<int, char> > > threadedCriticalTypes(threadNumber_);

#ifdef withOpenMP
  omp_lock_t writeLock;
  omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(int i = 0; i < (int) edgeList_->size(); i++){

    // avoid any processing if the abort signal is sent
    if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){

      int threadId = 0;
#ifdef withOpenMP
      threadId = omp_get_thread_num();
#endif
      
      // processing here!
      int pivotVertexId = (*edgeList_)[i].first;
      int otherExtremityId = (*edgeList_)[i].second;
     
      // A) compute the distance field
      double projectedPivotVertex[2];
      projectedPivotVertex[0] = uField[pivotVertexId];
      projectedPivotVertex[1] = vField[pivotVertexId];
      
      double projectedOtherVertex[2];
      projectedOtherVertex[0] = uField[otherExtremityId];
      projectedOtherVertex[1] = vField[otherExtremityId];
      
      double rangeEdge[2];
      rangeEdge[0] = projectedOtherVertex[0] - projectedPivotVertex[0];
      rangeEdge[1] = projectedOtherVertex[1] - projectedPivotVertex[1];
    
      double rangeNormal[2];
      rangeNormal[0] = -rangeEdge[1];
      rangeNormal[1] = rangeEdge[0];
      
      for(int j = 0; j < (int) (*edgeFans_)[i].size()/4; j++){
        for(int k = 0; k < 3; k++){
          
          int vertexId = (*edgeFans_)[i][j*4 + 1 + k];
          
          // we can compute the distance field (in the rage)
          double projectedVertex[2];
          projectedVertex[0] = uField[vertexId];
          projectedVertex[1] = vField[vertexId];
        
          double vertexRangeEdge[2];
          vertexRangeEdge[0] = projectedVertex[0] - projectedPivotVertex[0];
          vertexRangeEdge[1] = projectedVertex[1] - projectedPivotVertex[1];
          
          // signed distance: linear function of the dot product
          threadedDistanceField[threadId][vertexId] = 
            vertexRangeEdge[0]*rangeNormal[0] 
              + vertexRangeEdge[1]*rangeNormal[1];
        }
      }
      
      // B) compute critical points
      // watch out between local and global Ids
      // what I could do is to translate the ids from global to local
      // also, lots of things in there can be done out of the loop
      
      // in the loop
      char type = 
        threadedCriticalPoints[threadId].getCriticalType(pivotVertexId,
          (*edgeFanLinkEdgeLists_)[i]);
        
      if(type != -2){
        // -2: regular vertex
        threadedCriticalTypes[threadId].push_back(pair<int, char>(i, type));
      }
      
      // update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > advancedInfoMsg){
#ifdef withOpenMP
        omp_set_lock(&writeLock);
#endif
        if((wrapper_)
          &&(!(count % ((vertexNumber_)/10)))){
          wrapper_->updateProgress((count + 1.0)
            /vertexNumber_);
        }

        count++;
#ifdef withOpenMP
        omp_unset_lock(&writeLock);
#endif
      }
    }
  }
  
  // now merge the threaded lists
  for(int i = 0; i < threadNumber_; i++){
    for(int j = 0; j < (int) threadedCriticalTypes[i].size(); j++){
      jacobiSet.push_back(threadedCriticalTypes[i][j]);
    }
  }
  
#ifdef withOpenMP
  omp_destroy_lock(&writeLock);
#endif
 
  if(debugLevel_ >= Debug::infoMsg){
    int minimumNumber = 0, saddleNumber = 0, maximumNumber = 0, 
      monkeySaddleNumber = 0;
      
    for(int i = 0; i < (int) jacobiSet.size(); i++){
      switch(jacobiSet[i].second){
        case 0:
          minimumNumber++;
          break;
        case 1:
          saddleNumber++;
          break;
        case 2:
          maximumNumber++;
          break;
        case -1:
          monkeySaddleNumber++;
          break;
      }
    }
    
    {
      stringstream msg;
      msg << "[JacobiSet] Minimum edges: " << minimumNumber << endl;
      msg << "[JacobiSet] Saddle edges: " << saddleNumber << endl;
      msg << "[JacobiSet] Maximum edges: " << maximumNumber << endl;
      msg << "[JacobiSet] Multi-saddle edges: " << monkeySaddleNumber << endl;
      dMsg(cout, msg.str(), Debug::infoMsg);
    }
  }
  
  {
    stringstream msg;
    msg << "[JacobiSet] Data-set (" << edgeList_->size()
      << " edges) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    msg << "[JacobiSet] Jacobi edge rate: "
      << 100*(jacobiSet.size()/((double) edgeList_->size()))
      << "%" << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
  
  return 0;
}

template <class dataTypeU, class dataTypeV> 
  int JacobiSet<dataTypeU, dataTypeV>::perturbate(
    const dataTypeU &uEpsilon, const dataTypeV &vEpsilon) const{
      
#ifndef withKamikaze
  if(!uField_)
    return -1;
  if(!vField_)
    return -2;
  if(!vertexNumber_)
    return -3;
#endif
  
  dataTypeU *uField = (dataTypeU *) uField_;
  dataTypeV *vField = (dataTypeV *) vField_;
  
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < vertexNumber_; i++){
    // simulation of simplicity in 2 dimensions, need to use degree 2 polynoms
    
    uField[i] += i*uEpsilon;
    vField[i] += (i*vEpsilon)*(i*vEpsilon);
  }
      
  return 0;
}