#include                  <ScalarFieldCriticalPoints.h>

template <class dataType> 
  ScalarFieldCriticalPoints<dataType>::ScalarFieldCriticalPoints(){

  
  dimension_ = 0;
  vertexNumber_ = 0;
  scalarValues_ = NULL;
  vertexLinkEdgeLists_ = NULL;
  criticalPoints_ = NULL;
  sosOffsets_ = NULL;
 
//   threadNumber_ = 1;
}

template <class dataType> 
  ScalarFieldCriticalPoints<dataType>::~ScalarFieldCriticalPoints(){
  
}

template <class dataType> int ScalarFieldCriticalPoints<dataType>::execute(){

  // check the consistency of the variables -- to adapt
#ifndef withKamikaze
  if(!dimension_)
    return -1;
  if(!vertexNumber_)
    return -2;
  if(!scalarValues_)
    return -3;
  if(!vertexLinkEdgeLists_)
    return -4;
  if(!criticalPoints_)
    return -5;
#endif

  if(!sosOffsets_){
    // let's use our own local copy
    sosOffsets_ = &localSosOffSets_;
  }
  if((int) sosOffsets_->size() != vertexNumber_){
    Timer preProcess;
    sosOffsets_->resize(vertexNumber_);
    for(int i = 0; i < vertexNumber_; i++)
      (*sosOffsets_)[i] = i;
      
    {
      stringstream msg;
      msg << 
        "[ScalarFieldCriticalPoints] Offset pre-processing done in "
        << preProcess.getElapsedTime() << " s. Go!" << endl;
      dMsg(cout, msg.str(), timeMsg);
    }
  }
  
  Timer t;
  
  int count = 0;
 
  vector<char> vertexTypes(vertexNumber_);
  
#ifdef withOpenMP
  omp_lock_t writeLock;
  omp_init_lock(&writeLock);
#pragma omp parallel for num_threads(threadNumber_) 
#endif
  for(int i = 0; i < (int) vertexNumber_; i++){
  
    // avoid any processing if the abort signal is sent
    if((!wrapper_)||((wrapper_)&&(!wrapper_->needsToAbort()))){

      vertexTypes[i] = getCriticalType(i);
      
      // update the progress bar of the wrapping code -- to adapt
      if(debugLevel_ > advancedInfoMsg){
        if((wrapper_)
          &&(!(count % ((vertexNumber_)/10)))){
          wrapper_->updateProgress((count + 1.0)
            /vertexNumber_);
        }
#ifdef withOpenMP
        omp_set_lock(&writeLock);
#endif
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
  
  int minimumNumber = 0, maximumNumber = 0, saddleNumber = 0,
    oneSaddleNumber = 0, twoSaddleNumber = 0, monkeySaddleNumber = 0;
 
  // debug msg
  if(debugLevel_ >= Debug::infoMsg){
    if(dimension_ == 3){
      for(int i = 0; i < vertexNumber_; i++){
        switch(vertexTypes[i]){
          
          case 0:
            minimumNumber++;
            break;
            
          case 1:
            oneSaddleNumber++;
            break;
            
          case 2:
            twoSaddleNumber++;
            break;
          
          case 3:
            maximumNumber++;
            break;
            
          case -1:
            monkeySaddleNumber++;
            break;
        }
      }
    }
    else if(dimension_ == 2){
      for(int i = 0; i < vertexNumber_; i++){
        switch(vertexTypes[i]){
          
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
    }
    
    {
      stringstream msg;
      msg << "[ScalarFieldCriticalPoints] " << minimumNumber << " minima."
        << endl;
      if(dimension_ == 3){
        msg << "[ScalarFieldCriticalPoints] " << oneSaddleNumber
          << " 1-saddle(s)." << endl;
        msg << "[ScalarFieldCriticalPoints] " << twoSaddleNumber
          << " 2-saddle(s)." << endl;
      }
      if(dimension_ == 2){
        msg << "[ScalarFieldCriticalPoints] " << saddleNumber
          << " saddle(s)." << endl;
      }
      msg << "[ScalarFieldCriticalPoints] " << monkeySaddleNumber
        << " multi-saddle(s)." << endl;
      msg << "[ScalarFieldCriticalPoints] " << maximumNumber
        << " maxima." << endl;
        
      msg << "[ScalarFieldCriticalPoints] Euler characteristic approximation:";
//       if(monkeySaddleNumber){
//         msg << " approximation";
//       }
//       msg << ": ";
      if(dimension_ == 3){
        msg 
          << minimumNumber - oneSaddleNumber + twoSaddleNumber - maximumNumber;
      }
      if(dimension_ == 2){
        msg << minimumNumber - saddleNumber + maximumNumber;
      }
      msg << endl;
      dMsg(cout, msg.str(), Debug::infoMsg);
    }
  }
  
  // prepare the output
  criticalPoints_->clear();
  for(int i = 0; i < vertexNumber_; i++){
    if(vertexTypes[i] != -2){
      criticalPoints_->push_back(pair<int, char>(i, vertexTypes[i]));
    }
  }
  
  {
    stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Data-set (" << vertexNumber_
      << " vertices) processed in "
      << t.getElapsedTime() << " s. (" << threadNumber_
      << " thread(s))."
      << endl;
    dMsg(cout, msg.str(), 2);
  }
  
  return 0;
}

template <class dataType> char ScalarFieldCriticalPoints<dataType>
  ::getCriticalType(const int &vertexId, 
    const vector<pair<int, int> > &vertexLink) const{

  map<int, int> global2LowerLink, global2UpperLink;
  map<int, int>::iterator neighborIt;
  
  int lowerCount = 0, upperCount = 0;
  
  for(int i = 0; i < (int) vertexLink.size(); i++){
    
    int neighborId = vertexLink[i].first;
   
    // first vertex
    // lower link search
    if(isSosLowerThan(
      (*sosOffsets_)[neighborId], scalarValues_[neighborId],
      (*sosOffsets_)[vertexId], scalarValues_[vertexId])){
      
      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()){
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }
    
    // upper link
    if(isSosHigherThan(
      (*sosOffsets_)[neighborId], scalarValues_[neighborId],
      (*sosOffsets_)[vertexId], scalarValues_[vertexId])){
      
      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()){
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }
    
    // second vertex
    neighborId = vertexLink[i].second;
    
    // lower link search
    if(isSosLowerThan(
      (*sosOffsets_)[neighborId], scalarValues_[neighborId],
      (*sosOffsets_)[vertexId], scalarValues_[vertexId])){
      
      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()){
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }
    
    // upper link
    if(isSosHigherThan(
      (*sosOffsets_)[neighborId], scalarValues_[neighborId],
      (*sosOffsets_)[vertexId], scalarValues_[vertexId])){
      
      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()){
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }
  }
  
  if(debugLevel_ >= Debug::advancedInfoMsg){
    stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId 
      << " lower link (" << lowerCount << " vertices)" << endl;;
    
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId 
      << " upper link (" << upperCount << " vertices)" << endl;
    
    dMsg(cout, msg.str(), Debug::advancedInfoMsg);
  }
  
  if(!lowerCount){
    // minimum
    return 0;
  }
  if(!upperCount){
    // maximum
    return dimension_;
  }
  
  // so far 40% of the computation, that's ok.
 
  // now enumerate the connected components of the lower and upper links
  // NOTE: a breadth first search might be faster than a UF
  // if so, one would need the one-skeleton data structure, not the edge list
  vector<UnionFind> lowerSeeds(lowerCount);
  vector<UnionFind> upperSeeds(upperCount);
  vector<UnionFind *> lowerList(lowerCount);
  vector<UnionFind *> upperList(upperCount);
  for(int i = 0; i < (int) lowerList.size(); i++)
    lowerList[i] = &(lowerSeeds[i]);
  for(int i = 0; i < (int) upperList.size(); i++)
    upperList[i] = &(upperSeeds[i]);
  
  for(int i = 0; i < (int) vertexLink.size(); i++){
    
    int neighborId0 = vertexLink[i].first;
    int neighborId1 = vertexLink[i].second;
    
    // process the lower link
    if((isSosLowerThan((*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
      (*sosOffsets_)[vertexId], scalarValues_[vertexId]))
      &&
      (isSosLowerThan((*sosOffsets_)[neighborId1], scalarValues_[neighborId1],
      (*sosOffsets_)[vertexId], scalarValues_[vertexId]))){
      
      // both vertices are lower, let's add that edge and update the UF
      map<int, int>::iterator n0It = global2LowerLink.find(neighborId0);
      map<int, int>::iterator n1It = global2LowerLink.find(neighborId1);
    
      lowerList[n0It->second] = 
        makeUnion(lowerList[n0It->second], lowerList[n1It->second]);
      lowerList[n1It->second] = lowerList[n0It->second];
    }
    
    // process the upper link
    if((isSosHigherThan((*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
      (*sosOffsets_)[vertexId], scalarValues_[vertexId]))
      &&
      (isSosHigherThan((*sosOffsets_)[neighborId1], scalarValues_[neighborId1],
      (*sosOffsets_)[vertexId], scalarValues_[vertexId]))){
      
      // both vertices are lower, let's add that edge and update the UF
      map<int, int>::iterator n0It = global2UpperLink.find(neighborId0);
      map<int, int>::iterator n1It = global2UpperLink.find(neighborId1);
    
      upperList[n0It->second] = 
        makeUnion(upperList[n0It->second], upperList[n1It->second]);
      upperList[n1It->second] = upperList[n0It->second];
    }
  }
  
  // let's remove duplicates
  vector<UnionFind *>::iterator it;
  // update the UFs if necessary
  for(int i = 0; i < (int) lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(int i = 0; i < (int) upperList.size(); i++)
    upperList[i] = upperList[i]->find();
  
  sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));
  
  sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));
  
  if(debugLevel_ >= Debug::advancedInfoMsg){
    stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId
      << ": lowerLink-#CC=" << lowerList.size() 
      << " upperLink-#CC=" << upperList.size() << endl;
      
    dMsg(cout, msg.str(), Debug::advancedInfoMsg);
  }
  
  if((lowerList.size() == 1)&&(upperList.size() == 1))
    // regular point
    return -2;
  else{
    // saddles
    if(dimension_ == 2){
      if((lowerList.size() > 2)||(upperList.size() > 2)){
        // monkey saddle
        return -1;
      }
      else{
        // regular saddle
        return 1;
        // NOTE: you may have multi-saddles on the boundary in that 
        // configuration
        // to make this computation 100% correct, one would need to disambiguate
        // boundary from interior vertices
      }
    }
    else if(dimension_ == 3){
      if((lowerList.size() == 2)&&(upperList.size() == 1)){
        return 1;
      }
      else if((lowerList.size() == 1)&&(upperList.size() == 2)){
        return 2;
      }
      else{
        // monkey saddle
        return -1;
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }
  
  // -2: regular points
  return -2;
}