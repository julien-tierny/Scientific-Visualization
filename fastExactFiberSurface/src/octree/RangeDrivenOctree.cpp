/* 
 * file:                  RangeDrivenOctree.cpp
 * description:           RangeDrivenOctree processing package.
 * author:                Julien Tierny <julien.tierny@lip6.fr>.
 * date:                  March 2015.
 */

// TODO
// this class has been implemented recursively. this is bad programming.
// a non-recursive version will be implemented asap.

#include                  <RangeDrivenOctree.h>

template <class dataType0, class dataType1> 
  RangeDrivenOctree<dataType0, dataType1>::RangeDrivenOctree(){

  cellList_ = NULL;
  cellNumber_ = 0;
  
  leafMinimumCellNumber_ = 6;
  leafMinimumDomainVolumeRatio_ = 0.01;
  leafMinimumRangeAreaRatio_ = 0.01;
  
  pointList_ = NULL;
  
  u_ = NULL;
  v_ = NULL;
  
  vertexNumber_ = 0;
  rootId_ = 0; 
}

template <class dataType0, class dataType1> 
  RangeDrivenOctree<dataType0, dataType1>::~RangeDrivenOctree(){
  
}

template <class dataType0, class dataType1> 
  int RangeDrivenOctree<dataType0, dataType1>::build(){
  
  Timer t;
  Memory m;
  
  cellDomainBox_.resize(cellNumber_, vector<pair<float, float> >(3));
  
  cellRangeBox_.resize(cellNumber_);

  // WARNING: assuming tets only here
    
#ifdef withOpenMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < cellNumber_; i++){
    
    for(int j = 0; j < 3; j++){
      cellDomainBox_[i][j].first = FLT_MAX;
      cellDomainBox_[i][j].second = -FLT_MAX;
    }
    
    const long long int *cell = &(cellList_[5*i + 1]);
    
    for(int j = 0; j < 4; j++){
      
      // update the domain bounding box for that tet
      
      // update x-range
      if(pointList_[3*cell[j]] < cellDomainBox_[i][0].first)
        cellDomainBox_[i][0].first = pointList_[3*cell[j]];
      if(pointList_[3*cell[j]] > cellDomainBox_[i][0].second)
        cellDomainBox_[i][0].second = pointList_[3*cell[j]];
        
      // update y-range
      if(pointList_[3*cell[j] + 1] < cellDomainBox_[i][1].first)
        cellDomainBox_[i][1].first = pointList_[3*cell[j] + 1];
      if(pointList_[3*cell[j] + 1] > cellDomainBox_[i][1].second)
        cellDomainBox_[i][1].second = pointList_[3*cell[j] + 1];
      
      // update z-range
      if(pointList_[3*cell[j] + 2] < cellDomainBox_[i][2].first)
        cellDomainBox_[i][2].first = pointList_[3*cell[j] + 2];
      if(pointList_[3*cell[j] + 2] > cellDomainBox_[i][2].second)
        cellDomainBox_[i][2].second = pointList_[3*cell[j] + 2];
        
      // update the range bounding box
      if(!j){
        cellRangeBox_[i].first.first = u_[cell[j]];
        cellRangeBox_[i].first.second = u_[cell[j]];
        
        cellRangeBox_[i].second.first = v_[cell[j]];
        cellRangeBox_[i].second.second = v_[cell[j]];
      }
      else{
        
        // update u
        if(u_[cell[j]] < cellRangeBox_[i].first.first)
          cellRangeBox_[i].first.first = u_[cell[j]];
        if(u_[cell[j]] > cellRangeBox_[i].first.second)
          cellRangeBox_[i].first.second = u_[cell[j]];
        
        // update v
        if(v_[cell[j]] < cellRangeBox_[i].second.first)
          cellRangeBox_[i].second.first = v_[cell[j]];
        if(v_[cell[j]] > cellRangeBox_[i].second.second)
          cellRangeBox_[i].second.second = v_[cell[j]];
      }
    }
  }

  vector<int> domain(cellNumber_);
  for(int i = 0; i < cellNumber_; i++)
    domain[i] = i;
  
  // get global bBoxes
  vector<pair<float, float> > domainBox(3);
  pair<pair<dataType0, dataType0>, pair<dataType1, dataType1> > rangeBox;
  
  for(int i = 0; i < vertexNumber_; i++){
    
    // domain one
    for(int j = 0; j < 3; j++){
      if(!i){
        domainBox[j].first = domainBox[j].second = pointList_[3*i + j];
      }
      else{
        if(pointList_[3*i + j] < domainBox[j].first)
          domainBox[j].first = pointList_[3*i + j];
        if(pointList_[3*i + j] > domainBox[j].second)
          domainBox[j].second = pointList_[3*i + j];
      }
    }
    
    if(!i){
      rangeBox.first.first = rangeBox.first.second = u_[i];
      rangeBox.second.first = rangeBox.second.second = v_[i];
    }
    else{
      if(u_[i] < rangeBox.first.first)
        rangeBox.first.first = u_[i];
      if(u_[i] > rangeBox.first.second)
        rangeBox.first.second = u_[i];
      
      if(v_[i] < rangeBox.second.first)
        rangeBox.second.first = v_[i];
      if(v_[i] > rangeBox.second.second)
        rangeBox.second.second = v_[i];
    }
  }
  
  rangeArea_ = (rangeBox.first.second - rangeBox.first.first)
    *(rangeBox.second.second - rangeBox.second.first);
  domainVolume_ = (domainBox[0].second - domainBox[0].first)
    *(domainBox[1].second - domainBox[1].first)
    *(domainBox[2].second - domainBox[2].first);
  
  // special case for tets obtained from regular grid subdivision (assuming 6)
  if(leafMinimumCellNumber_ < 6)
    leafMinimumCellNumber_ = 6;
  
  leafMinimumDomainVolumeRatio_ = (1.0/((float)cellNumber_))/2.0;
  
  {
    stringstream msg;
    msg << "[RangeDrivenOctree] Range area ratio: "
      << leafMinimumRangeAreaRatio_
      << endl;
    dMsg(cout, msg.str(), 4);
  }
  buildNode(domain, domainBox, rangeBox, rootId_);
 
  {
    stringstream msg;
    msg << "[RangeDrivenOctree] Octree built in "
      << t.getElapsedTime() << " s." << endl;
    dMsg(cout, msg.str(), 2);
  }
  {
    stringstream msg;
    msg << "[RangeDrivenOctree] Memory: "
      << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), 2);
  }
  
  // debug
  stats(cout);
  // end of debug
  
  return 0;    
}

template <class dataType0, class dataType1> 
  int RangeDrivenOctree<dataType0, dataType1>::buildNode(
    const vector<int> &cellList,
    const vector<pair<float, float> > &domainBox,
    const pair<pair<dataType0, dataType0>, pair<dataType1, dataType1> > 
      &rangeBox, 
    int &nodeId){

  nodeId = nodeList_.size();
  
  nodeList_.resize(nodeList_.size() + 1);

  nodeList_.back().rangeBox_ = rangeBox;
  nodeList_.back().domainBox_ = domainBox;
  
  float rangeArea = (rangeBox.first.second - rangeBox.first.first)
    *(rangeBox.second.second - rangeBox.second.first);
    
  float domainVolume = (domainBox[0].second - domainBox[0].first)
    *(domainBox[1].second - domainBox[1].first)
    *(domainBox[2].second - domainBox[2].first);
  
  if(((int) cellList.size() > leafMinimumCellNumber_)
    &&(rangeArea > leafMinimumRangeAreaRatio_*rangeArea_)
    &&(domainVolume > leafMinimumDomainVolumeRatio_*domainVolume_)){
    
    nodeList_.back().childList_.resize(8);
  
    vector<vector<pair<float, float> > > childDomainBox(8);
    for(int i = 0; i < (int) childDomainBox.size(); i++){
      childDomainBox[i].resize(3);
    }
    vector<pair<pair<dataType0, dataType0>, pair<dataType1, dataType1> > > 
      childRangeBox(8);
    vector<vector<int> > childCellList(8);

    float midX = domainBox[0].first 
      + (domainBox[0].second - domainBox[0].first)/2.0;
    float midY = domainBox[1].first
      + (domainBox[1].second - domainBox[1].first)/2.0;
    float midZ = domainBox[2].first
      + (domainBox[2].second - domainBox[2].first)/2.0;
     
    // 0 - - -
    childDomainBox[0][0].first = domainBox[0].first;
    childDomainBox[0][0].second = midX;
    childDomainBox[0][1].first = domainBox[1].first;
    childDomainBox[0][1].second = midY;
    childDomainBox[0][2].first = domainBox[2].first;
    childDomainBox[0][2].second = midZ;
   
    // 1 - - +
    childDomainBox[1][0].first = domainBox[0].first;
    childDomainBox[1][0].second = midX;
    childDomainBox[1][1].first = domainBox[1].first;
    childDomainBox[1][1].second = midY;
    childDomainBox[1][2].first = midZ;
    childDomainBox[1][2].second = domainBox[2].second;
    
    // 2 - + -
    childDomainBox[2][0].first = domainBox[0].first;
    childDomainBox[2][0].second = midX;
    childDomainBox[2][1].first = midY;
    childDomainBox[2][1].second = domainBox[1].second;
    childDomainBox[2][2].first = domainBox[2].first;
    childDomainBox[2][2].second = midZ;
    
    // 3 - + + 
    childDomainBox[3][0].first = domainBox[0].first;
    childDomainBox[3][0].second = midX;
    childDomainBox[3][1].first = midY;
    childDomainBox[3][1].second = domainBox[1].second;
    childDomainBox[3][2].first = midZ;
    childDomainBox[3][2].second = domainBox[2].second;
    
    // 4 + - -
    childDomainBox[4][0].first = midX;
    childDomainBox[4][0].second = domainBox[0].second;
    childDomainBox[4][1].first = domainBox[1].first;
    childDomainBox[4][1].second = midY;
    childDomainBox[4][2].first = domainBox[2].first;
    childDomainBox[4][2].second = midZ;
    
    // 5 + - +
    childDomainBox[5][0].first = midX;
    childDomainBox[5][0].second = domainBox[0].second;
    childDomainBox[5][1].first = domainBox[1].first;
    childDomainBox[5][1].second = midY;
    childDomainBox[5][2].first = midZ;
    childDomainBox[5][2].second = domainBox[2].second;
    
    // 6 + + -
    childDomainBox[6][0].first = midX;
    childDomainBox[6][0].second = domainBox[0].second;
    childDomainBox[6][1].first = midY;
    childDomainBox[6][1].second = domainBox[1].second;
    childDomainBox[6][2].first = domainBox[2].first;
    childDomainBox[6][2].second = midZ;
    
    // 7 + + +
    childDomainBox[7][0].first = midX;
    childDomainBox[7][0].second = domainBox[0].second;
    childDomainBox[7][1].first = midY;
    childDomainBox[7][1].second = domainBox[1].second;
    childDomainBox[7][2].first = midZ;
    childDomainBox[7][2].second = domainBox[2].second;
    
    for(int i = 0; i < (int) cellList.size(); i++){
      
      int childId = 0;
      
      for(int j = 0; j < 8; j++){
        if((cellDomainBox_[cellList[i]][0].first 
          >= childDomainBox[j][0].first)
          &&(cellDomainBox_[cellList[i]][0].first
          < childDomainBox[j][0].second)
          &&
          (cellDomainBox_[cellList[i]][1].first 
          >= childDomainBox[j][1].first)
          &&(cellDomainBox_[cellList[i]][1].first
          < childDomainBox[j][1].second)
          &&
          (cellDomainBox_[cellList[i]][2].first 
          >= childDomainBox[j][2].first)
          &&(cellDomainBox_[cellList[i]][2].first
          < childDomainBox[j][2].second)){
        
          childId = j;
          break;
        }
      }
      
      // update child's range box
      if(childCellList[childId].empty()){
        childRangeBox[childId].first.first 
          = cellRangeBox_[cellList[i]].first.first;
        childRangeBox[childId].first.second
          = cellRangeBox_[cellList[i]].first.second;
          
        childRangeBox[childId].second.first 
          = cellRangeBox_[cellList[i]].second.first;
        childRangeBox[childId].second.second
          = cellRangeBox_[cellList[i]].second.second;
      }
      else{
        if(cellRangeBox_[cellList[i]].first.first 
          < childRangeBox[childId].first.first){
          childRangeBox[childId].first.first 
            = cellRangeBox_[cellList[i]].first.first;
        }
        if(cellRangeBox_[cellList[i]].first.second
          > childRangeBox[childId].first.second){
          childRangeBox[childId].first.second 
            = cellRangeBox_[cellList[i]].first.second;
        }
        
        if(cellRangeBox_[cellList[i]].second.first 
          < childRangeBox[childId].second.first){
          childRangeBox[childId].second.first 
            = cellRangeBox_[cellList[i]].second.first;
        }
        if(cellRangeBox_[cellList[i]].second.second
          > childRangeBox[childId].second.second){
          childRangeBox[childId].second.second 
            = cellRangeBox_[cellList[i]].second.second;
        }
      }
    
      childCellList[childId].push_back(cellList[i]);
      
    }
    
    for(int i = 0; i < 8; i++){
      buildNode(childCellList[i], 
        childDomainBox[i], childRangeBox[i], nodeList_[nodeId].childList_[i]);
    }
  }
  else{
    // leaf
    nodeList_[nodeId].cellList_ = cellList;
  }
  
  return 0;
}

template<class dataType0, class dataType1>
  inline int RangeDrivenOctree<dataType0, dataType1>::getTet2NodeMap(
    vector<int> &map, const bool &forSegmentation) const{

  vector<int> randomMap;
  if(forSegmentation){
    randomMap.resize(nodeList_.size());
    for(int i = 0; i < (int) randomMap.size(); i++){
      randomMap[i] = rand()%(randomMap.size());
    }
  }
      
  map.resize(cellNumber_);
  for(int i = 0; i < (int) nodeList_.size(); i++){
    for(int j = 0; j < (int) nodeList_[i].cellList_.size(); j++){
      if(forSegmentation){
       map[nodeList_[i].cellList_[j]] = randomMap[i]; 
      }
      else{
        map[nodeList_[i].cellList_[j]] = i;
      }
    }
  }
    
  return 0;
}


template<class dataType0, class dataType1>
  inline int RangeDrivenOctree<dataType0, dataType1>::rangeSegmentQuery(
    const double *p0, const double *p1,
    vector<int> &cellList) const{

  Timer t;
  
  queryResultNumber_ = 0;
  cellList.clear();
  
  int ret = rangeSegmentQuery(p0, p1, rootId_, cellList);
  
  {
    stringstream msg;
    msg << "[RangeDrivenOctree] Query done in "
      << t.getElapsedTime() << " s. ("
      << queryResultNumber_ << " non-empty leaves, "
      << cellList.size() << " cells)" << endl;
    dMsg(cout, msg.str(), 2);
  }
  
  return ret;
}

template<class dataType0, class dataType1>
  inline int RangeDrivenOctree<dataType0, dataType1>::rangeSegmentQuery(
    const double *p0, const double *p1,
    const int &nodeId,
    vector<int> &cellList) const{

  // check for intersection for each segment of the range bounding box 
  double q0[2], q1[2];
 
  // bottom range segment (min, min) (max, min)
  q0[0] = nodeList_[nodeId].rangeBox_.first.first;
  q0[1] = nodeList_[nodeId].rangeBox_.second.first;
  
  q1[0] = nodeList_[nodeId].rangeBox_.first.second;
  q1[1] = q0[1];
  
  if(segmentIntersection(p0, p1, q0, q1)){
    if(nodeList_[nodeId].childList_.size()){
      for(int i = 0; i < (int) nodeList_[nodeId].childList_.size(); i++){
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    }
    else{
      // terminal leaf
      // return our cells
      {
        stringstream msg;
        msg << "[RangeDrivenOctree] Node #"
          << nodeId << " returns its "
          << nodeList_[nodeId].cellList_.size() 
          << " cell(s)." << endl;
        dMsg(cout, msg.str(), 4);
      }
      cellList.insert(cellList.end(), 
        nodeList_[nodeId].cellList_.begin(), nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }
  
  // right segment (max, min) (max, max)
  q0[0] = nodeList_[nodeId].rangeBox_.first.second;
  q0[1] = nodeList_[nodeId].rangeBox_.second.first;
  
  q1[0] = nodeList_[nodeId].rangeBox_.first.second;
  q1[1] = nodeList_[nodeId].rangeBox_.second.second;
  
  if(segmentIntersection(p0, p1, q0, q1)){
    if(nodeList_[nodeId].childList_.size()){
      for(int i = 0; i < (int) nodeList_[nodeId].childList_.size(); i++){
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    }
    else{
      // terminal leaf
      // return our cells
      {
        stringstream msg;
        msg << "[RangeDrivenOctree] Node #"
          << nodeId << " returns its "
          << nodeList_[nodeId].cellList_.size() 
          << " cell(s)." << endl;
        dMsg(cout, msg.str(), 4);
      }
      cellList.insert(cellList.end(), 
        nodeList_[nodeId].cellList_.begin(), nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }
  
  // top segment (min, max) (max, max)
  q0[0] = nodeList_[nodeId].rangeBox_.first.first;
  q0[1] = nodeList_[nodeId].rangeBox_.second.second;
  
  q1[0] = nodeList_[nodeId].rangeBox_.first.second;
  q1[1] = nodeList_[nodeId].rangeBox_.second.second;
  
  if(segmentIntersection(p0, p1, q0, q1)){
    if(nodeList_[nodeId].childList_.size()){
      for(int i = 0; i < (int) nodeList_[nodeId].childList_.size(); i++){
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    }
    else{
      // terminal leaf
      // return our cells
      {
        stringstream msg;
        msg << "[RangeDrivenOctree] Node #"
          << nodeId << " returns its "
          << nodeList_[nodeId].cellList_.size() 
          << " cell(s)." << endl;
        dMsg(cout, msg.str(), 4);
      }
      cellList.insert(cellList.end(), 
        nodeList_[nodeId].cellList_.begin(), nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }
  
  // left segment (min, min) (min, max)
  q0[0] = nodeList_[nodeId].rangeBox_.first.first;
  q0[1] = nodeList_[nodeId].rangeBox_.second.first;
  
  q1[0] = nodeList_[nodeId].rangeBox_.first.first;
  q1[1] = nodeList_[nodeId].rangeBox_.second.second;
  
  if(segmentIntersection(p0, p1, q0, q1)){
    if(nodeList_[nodeId].childList_.size()){
      for(int i = 0; i < (int) nodeList_[nodeId].childList_.size(); i++){
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    }
    else{
      // terminal leaf
      // return our cells
      {
        stringstream msg;
        msg << "[RangeDrivenOctree] Node #"
          << nodeId << " returns its "
          << nodeList_[nodeId].cellList_.size() 
          << " cell(s)." << endl;
        dMsg(cout, msg.str(), 4);
      }
      cellList.insert(cellList.end(), 
        nodeList_[nodeId].cellList_.begin(), nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }
  
  // is the segment completely included in the range bounding box?
  if((p0[0] >= nodeList_[nodeId].rangeBox_.first.first)
    &&(p0[0] < nodeList_[nodeId].rangeBox_.first.second)
    &&(p0[1] >= nodeList_[nodeId].rangeBox_.second.first)
    &&(p0[1] < nodeList_[nodeId].rangeBox_.second.second)){
    
    // p0 is in there
    if(nodeList_[nodeId].childList_.size()){
      for(int i = 0; i < (int) nodeList_[nodeId].childList_.size(); i++){
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    }
    else{
      // terminal leaf
      // return our cells
      {
        stringstream msg;
        msg << "[RangeDrivenOctree] Node #"
          << nodeId << " returns its "
          << nodeList_[nodeId].cellList_.size() 
          << " cell(s)." << endl;
        dMsg(cout, msg.str(), 4);
      }
      cellList.insert(cellList.end(), 
        nodeList_[nodeId].cellList_.begin(), nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }
  if((p1[0] >= nodeList_[nodeId].rangeBox_.first.first)
    &&(p1[0] < nodeList_[nodeId].rangeBox_.first.second)
    &&(p1[1] >= nodeList_[nodeId].rangeBox_.second.first)
    &&(p1[1] < nodeList_[nodeId].rangeBox_.second.second)){
    
    // p1 is in there
    if(nodeList_[nodeId].childList_.size()){
      for(int i = 0; i < (int) nodeList_[nodeId].childList_.size(); i++){
        rangeSegmentQuery(p0, p1, nodeList_[nodeId].childList_[i], cellList);
      }
      return 0;
    }
    else{
      // terminal leaf
      // return our cells
      {
        stringstream msg;
        msg << "[RangeDrivenOctree] Node #"
          << nodeId << " returns its "
          << nodeList_[nodeId].cellList_.size() 
          << " cell(s)." << endl;
        dMsg(cout, msg.str(), 4);
      }
      cellList.insert(cellList.end(), 
        nodeList_[nodeId].cellList_.begin(), nodeList_[nodeId].cellList_.end());
      queryResultNumber_++;
      return 0;
    }
  }
    
  return 0;
}

template<class dataType0, class dataType1>
  int RangeDrivenOctree<dataType0, dataType1>::statNode(const int &nodeId,
    ostream &stream){

  stream << "[RangeDrivenOctree]" << endl;
  stream << "[RangeDrivenOctree] Node #" << nodeId << endl;
  stream << "[RangeDrivenOctree]   Domain box: ["
    << nodeList_[nodeId].domainBox_[0].first << " "
    << nodeList_[nodeId].domainBox_[0].second << "] ["
    << nodeList_[nodeId].domainBox_[1].first << " "
    << nodeList_[nodeId].domainBox_[1].second << "] ["
    << nodeList_[nodeId].domainBox_[2].first << " "
    << nodeList_[nodeId].domainBox_[2].second << "] "
    << " volume=" << 
      (nodeList_[nodeId].domainBox_[0].second 
        - nodeList_[nodeId].domainBox_[0].first)
      *(nodeList_[nodeId].domainBox_[1].second 
        - nodeList_[nodeId].domainBox_[1].first)
      *(nodeList_[nodeId].domainBox_[2].second 
        - nodeList_[nodeId].domainBox_[2].first)
    << " threshold=" << leafMinimumDomainVolumeRatio_*domainVolume_
    << endl;
  stream << "[RangeDrivenOctree]   Range box: ["
    << nodeList_[nodeId].rangeBox_.first.first << " "
    << nodeList_[nodeId].rangeBox_.first.second << "] ["
    << nodeList_[nodeId].rangeBox_.second.first << " "
    << nodeList_[nodeId].rangeBox_.second.second << "] "
    << " area=" <<
      (nodeList_[nodeId].rangeBox_.first.second
        - nodeList_[nodeId].rangeBox_.first.first)
      *(nodeList_[nodeId].rangeBox_.second.second
        - nodeList_[nodeId].rangeBox_.second.first)
    << " threshold=" << leafMinimumRangeAreaRatio_*rangeArea_
    << endl;
  stream << "[RangeDrivenOctree] Number of cells: " 
    << nodeList_[nodeId].cellList_.size() << endl;
    
  return 0;
}

template<class dataType0, class dataType1> 
  int RangeDrivenOctree<dataType0, dataType1>::stats(ostream &stream){
  
  int leafNumber = 0, nonEmptyLeafNumber = 0,
      minCellNumber = -1, maxCellNumber = -1, storedCellNumber = 0;
  float averageCellNumber = 0;
  int maxCellId = 0;

  for(int i = 0; i < (int) nodeList_.size(); i++){
    if(!nodeList_[i].childList_.size()){
      // leaf
      leafNumber++;
      if(nodeList_[i].cellList_.size()){
        nonEmptyLeafNumber++;
        storedCellNumber += nodeList_[i].cellList_.size();
        
        averageCellNumber += nodeList_[i].cellList_.size();
        if((minCellNumber == -1)
          ||(nodeList_[i].cellList_.size() < minCellNumber))
          minCellNumber = nodeList_[i].cellList_.size();
        if((maxCellNumber == -1)
          ||(nodeList_[i].cellList_.size() > maxCellNumber)){
          maxCellNumber = nodeList_[i].cellList_.size();
          maxCellId = i;
        }
      }
    }
  }
  averageCellNumber /= nonEmptyLeafNumber;
  
  stream << "[RangeDrivenOctree] Domain volume: " << domainVolume_ << endl;
  stream << "[RangeDrivenOctree] Range area: " << rangeArea_ << endl;
  stream << "[RangeDrivenOctree] Leaf number: " 
    << leafNumber << endl;
  stream << "[RangeDrivenOctree] Non empty leaf number: "
    << nonEmptyLeafNumber << endl;
  stream << "[RangeDrivenOctree] Average cell number: "
    << averageCellNumber << endl;
  stream << "[RangeDrivenOctree] Min cell number: "
    << minCellNumber << endl;
  stream << "[RangeDrivenOctree] Max cell number: "
    << maxCellNumber << endl;
  stream << "[RangeDrivenOctree] Stored cell number: "
    << storedCellNumber << " [input=" << cellNumber_ << "]" << endl;
  
  stream << "[RangeDrivenOctree] Max-cell nodeId: " << maxCellId << endl;
 
  if(debugLevel_ > 5){
    for(int i = 0; i < (int) nodeList_.size(); i++){
      if(nodeList_[i].cellList_.size())
        statNode(i, stream);
    }
  }

  return 0;  
}

template<class dataType0, class dataType1> 
  inline bool RangeDrivenOctree<dataType0, dataType1>::segmentIntersection(
    const double *p0, const double *p1, 
    const double *q0, const double *q1) const{

  // NOTE: [q0, q1] is axis aligned
  bool horizontal = true;
  double x = 0, y = 0;
  
  if(q0[0] == q1[0]) horizontal = false;
  
  double denP = p1[0] - p0[0];
  if(!denP) denP = DBL_EPSILON;
  
  double P = (p1[1] - p0[1])/denP;
  if(!P) P = DBL_EPSILON;
  double bP = p1[1] - P*p1[0];
  
  if(horizontal){
    y = q0[1];
    x = (y - bP)/P;
  }
  else{
    x = q0[0];
    y = P*x + bP;
  }
  
  // does the intersection lands on the range segment?
  if(horizontal){
    if((x < fmin(q0[0], q1[0]))||(x > fmax(q0[0], q1[0]))){
      return false;
    }
  }
  else{
     if((y < fmin(q0[1], q1[1]))||(y > fmax(q0[1], q1[1]))){
      return false;
    }
  }
  
  // does it lands on the polygon segment's projection?
  if((x < fmin(p0[0], p1[0])||(x > fmax(p0[0], p1[0])))
    &&((y < fmin(p0[1], p1[1])||(y > fmax(p0[1], p1[1]))))){
    return false;
  }
  
  return true;
}
