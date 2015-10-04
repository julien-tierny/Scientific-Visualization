/*
 * file:                TriangleMesh.cpp
 * description:         Lightweight class for PL 2-manifold handling.
 * author:              (C) Julien Tierny <tierny@telecom-paristech.fr>
 * date:                March 2012.
 *
 * This program is free software: you can redistribute it and/or modify it 
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version. 
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License 
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include                <TriangleMesh.h>

// Vertex implementation
Vertex::Vertex(){
  mesh_ = NULL;
  point_ = NULL;
}

Vertex::Vertex(const Vertex &other){

  if(&other != this){
    point_ = other.point_;
    edgeIds_ = other.edgeIds_;
    triangleIds_ = other.triangleIds_;

    // TODO: complete with the rest of the structure

    mesh_ = other.mesh_;
  }
}

Vertex::Vertex(const Vertex &other, TriangleMesh *mesh){

  if(&other != this){
    point_ = other.point_;
    edgeIds_ = other.edgeIds_;
    triangleIds_ = other.triangleIds_;
    mesh_ = mesh;
  }
}

Vertex::Vertex(const vector<double> &point){

  if(point.size() != 3)
    dMsg(cerr,
      "[TriangleMesh] Ignoring extra coordinates for non 3D point.\n", 10);
}

Vertex::~Vertex(){}

Vertex& Vertex::operator=(const Vertex &other){

  if(&other != this){
    if(point_)
      for(int i = 0; i < 3; i++)
        (*point_)[i] = (*other.point_)[i];
    edgeIds_ = other.edgeIds_;
    triangleIds_ = other.triangleIds_;

    mesh_ = other.mesh_;
  }

  return *this;
}

void Vertex::operator<<(istream &f){

  if(point_)
    for(int i = 0; i < 3; i++)
      f >> (*point_)[i];
}

void Vertex::operator>>(ostream &f) const{

  if(point_){
    for(int i = 0; i < 3; i++)
      f << " " << (*point_)[i];
    f << endl;
  }
}

int Vertex::exportToVRML(ostream &f) const{

  if(point_){
    f << "\t\t\t\t\t\t" 
      << (*point_)[0] << " "
      << (*point_)[1] << " "
      << (*point_)[2] << endl;
  }
  return 0;
}

const Edge* Vertex::getEdge(const int &vertexEdgeId) const{

  if((vertexEdgeId < 0)||(vertexEdgeId >= edgeIds_.size())){
    stringstream msg;
    msg << "[Vertex] Vertex edge Id #" << vertexEdgeId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return NULL;
  }

  if(mesh_)
    return mesh_->getEdge(edgeIds_[vertexEdgeId]);
  else 
    return NULL;
}

int Vertex::getEdgeId(const int &vertexEdgeId, int &edgeId) const{

  if((vertexEdgeId < 0)||(vertexEdgeId >= edgeIds_.size())){
    stringstream msg;
    msg << "[Vertex] Vertex edge Id #" << vertexEdgeId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return -1;
  }

  edgeId = edgeIds_[vertexEdgeId];

  return 0;
}

const Edge* Vertex::getLinkEdge(const int &linkEdgeId) const{

  if((linkEdgeId < 0)||(linkEdgeId >= linkEdgeIds_.size())){
    stringstream msg;
    msg << "[Vertex] Star edge Id #" << linkEdgeId << " out of range."  << endl;
    dMsg(cerr, msg.str(), 10);
    return NULL;
  }

  if(mesh_) return mesh_->getEdge(linkEdgeIds_[linkEdgeId]);
  else return NULL;
}

int Vertex::getLinkEdgeId(const int &linkEdgeId, int &edgeId) const{

  if((linkEdgeId < 0)||(linkEdgeId >= linkEdgeIds_.size())){
    stringstream msg;
    msg << "[Vertex] Star edge Id #" << linkEdgeId << " out of range."  << endl;
    dMsg(cerr, msg.str(), 10);
    return -1;
  }

  edgeId = linkEdgeIds_[linkEdgeId];

  return 0;
}

int Vertex::getTriangleId(const int &vertexTriangleId, int &triangleId) const{

  if((vertexTriangleId < 0)||(vertexTriangleId >= triangleIds_.size())){
    stringstream msg;
    msg << "[Vertex] Vertex triangle Id #" << vertexTriangleId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return -1;
  }

  triangleId = triangleIds_[vertexTriangleId];

  return 0;
}


// Edge implementation
Edge::Edge(){

  vertexIds_[0] = -1;
  vertexIds_[1] = -1;

  triangleIds_[0] = -1;
  triangleIds_[1] = -1;

  mesh_ = NULL;
}

Edge::Edge(const Edge &other){
  if(&other != this){
    for(int i = 0; i < 2; i++)
      vertexIds_[i] = other.vertexIds_[i];
    for(int i = 0; i < 2; i++)
      triangleIds_[i] = other.triangleIds_[i];

    mesh_ = other.mesh_;
  }
}

Edge::Edge(const Edge &other, TriangleMesh *mesh){

  if(&other != this){
    for(int i = 0; i < 2; i++)
      vertexIds_[i] = other.vertexIds_[i];
    for(int i = 0; i < 2; i++)
      triangleIds_[i] = other.triangleIds_[i];

    mesh_ = mesh;
  }

  mesh_ = mesh;
}

Edge::~Edge(){}
    
Edge& Edge::operator=(const Edge &other){

  if(&other != this){
    for(int i = 0; i < 2; i++)
      vertexIds_[i] = other.vertexIds_[i];
    for(int i = 0; i < 2; i++)
      triangleIds_[i] = other.triangleIds_[i];

    mesh_ = other.mesh_;
  }

  return *this;
}

void Edge::operator<<(istream &f){}

void Edge::operator>>(ostream &f) const{}

const Triangle* Edge::getTriangle(const int &edgeTriangleId) const{

   if((edgeTriangleId < 0)||(edgeTriangleId >= 2)){
    stringstream msg;
    msg << "[Edge] Edge vertex Id #" << edgeTriangleId << " out of range." 
      << endl;
    dMsg(cerr, msg.str(), 10);

    return NULL;
  }

  if(mesh_)
    return mesh_->getTriangle(triangleIds_[edgeTriangleId]);
  else 
    return NULL;
 
  return 0;
}

int Edge::getTriangleId(const int &edgeTriangleId, int &triangleId) const{

   if((edgeTriangleId < 0)||(edgeTriangleId >= 2)){
    stringstream msg;
    msg << "[Edge] Edge vertex Id #" << edgeTriangleId << " out of range." 
      << endl;
    dMsg(cerr, msg.str(), 10);

    return -1;
  }

  triangleId = triangleIds_[edgeTriangleId];
 
  return 0;
}

const Vertex* Edge::getVertex(const int &edgeVertexId) const{

  if((edgeVertexId < 0)||(edgeVertexId >= 2)){
    stringstream msg;
    msg << "[Edge] Edge vertex Id #" << edgeVertexId << " out of range." 
      << endl;
    dMsg(cerr, msg.str(), 10);

    return NULL;
  }

  if(mesh_)
    return mesh_->getVertex(vertexIds_[edgeVertexId]);
  else 
    return NULL;
}

int Edge::getVertexId(const int &edgeVertexId, int &vertexId) const{

  if((edgeVertexId < 0)||(edgeVertexId >= 2)){
    stringstream msg;
    msg << "[Edge] Edge vertex Id #" << edgeVertexId << " out of range." 
      << endl;
    dMsg(cerr, msg.str(), 10);

    return -1;
  }

  vertexId = vertexIds_[edgeVertexId];

  return 0;
}

// Triangle implementation
Triangle::Triangle(){

  vertexIds_ = NULL;
  mesh_ = NULL;
}

Triangle::Triangle(const Triangle &other){

  if(&other != this){
    vertexIds_ = other.vertexIds_;
    for(int i = 0; i < 3; i++)
      edgeIds_[i] = other.edgeIds_[i];
    for(int i = 0; i < 3; i++)
      triangleIds_[i] = other.triangleIds_[i];

    mesh_ = other.mesh_;
  }
}

Triangle::Triangle(const Triangle &other, TriangleMesh *mesh){

  if(&other != this){
    vertexIds_ = other.vertexIds_;
    for(int i = 0; i < 3; i++)
      edgeIds_[i] = other.edgeIds_[i];
    for(int i = 0; i < 3; i++)
      triangleIds_[i] = other.triangleIds_[i];

    mesh_ = mesh;
  }
}

Triangle::~Triangle(){
  
  for(int i = 0; i < 3; i++){
    if(vertexIds_)
      (*vertexIds_)[i] = -1;
    edgeIds_[i] = -1;
    triangleIds_[i] = -1;
  }

}

Triangle::Triangle(const vector<int> &vertexList){

  if(vertexList.size() != 3)
    dMsg(cerr,
      "[TriangleMesh] Ignoring extra vertices for non triangular face.\n", 10);

  if(vertexIds_)
    for(int i = 0; (i < vertexList.size())||(i < 3); i++)
      (*vertexIds_)[i] = vertexList[i];
}

Triangle& Triangle::operator=(const Triangle &other){
  
  if(&other != this){
    if(vertexIds_)
      for(int i = 0; i < 3; i++)
        (*vertexIds_)[i] = (*other.vertexIds_)[i];
    for(int i = 0; i < 3; i++)
      edgeIds_[i] = other.edgeIds_[i];
    for(int i = 0; i < 3; i++)
      triangleIds_[i] = other.triangleIds_[i];

    mesh_ = other.mesh_;
  }

  return *this;
}

void Triangle::operator<<(istream &f){
  
  int vertexNumber;

  f >> vertexNumber;

  if(vertexNumber != 3)
    dMsg(cerr, 
      "[TriangleMesh] Ignoring extra vertices for non triangular face.\n", 10);

  if(vertexIds_)
    for(int i = 0; (i < vertexNumber)||(i < 3); i++)
      f >> (*vertexIds_)[i];
}

void Triangle::operator>>(ostream &f) const{

  if(!vertexIds_) return;

  f << "3";
  for(int i = 0; i < 3; i++){
    f << " " << (*vertexIds_)[i];
  }
  f << endl;
}

int Triangle::exportToVRML(ostream &f) const{

  if(!vertexIds_) return -1;

  f << "\t\t\t\t\t"
    << (*vertexIds_)[0] << " "
    << (*vertexIds_)[1] << " "
    << (*vertexIds_)[2] << " -1" << endl;
  
  return 0;
}

const Edge* Triangle::getEdge(const int &triangleEdgeId) const{

  if((triangleEdgeId < 0)||(triangleEdgeId >= 3)){
    stringstream msg;
    msg << "[Triangle] Triangle edge Id #" << triangleEdgeId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);

    return NULL;
  }

  if(mesh_)
    return mesh_->getEdge(edgeIds_[triangleEdgeId]);
  else
    return NULL;

  return 0;
}

int Triangle::getEdgeId(const int &triangleEdgeId, int &edgeId) const{

  if((triangleEdgeId < 0)||(triangleEdgeId >= 3)){
    stringstream msg;
    msg << "[Triangle] Triangle edge Id #" << triangleEdgeId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);

    return -2;
  }

  edgeId = edgeIds_[triangleEdgeId];

  return 0;
}

const Triangle* Triangle::getTriangle(const int &triangleTriangleId) const{
   
  if((triangleTriangleId < 0)||(triangleTriangleId >= 3)){
    stringstream msg;
    msg << "[Triangle] Triangle triangle Id #" << triangleTriangleId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);

    return NULL;
  }

  if(mesh_)
    return mesh_->getTriangle(triangleIds_[triangleTriangleId]);
  else 
    return NULL;

  return 0;
}

int Triangle::getTriangleId(const int &triangleTriangleId, 
  int &triangleId) const{

  if((triangleTriangleId < 0)||(triangleTriangleId >= 3)){
    stringstream msg;
    msg << "[Triangle] Triangle triangle Id #" << triangleTriangleId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);

    return -2;
  }

  triangleId = triangleIds_[triangleTriangleId];

  return 0;
}

const Vertex* Triangle::getVertex(const int &triangleVertexId) const{

  if((triangleVertexId < 0)||(triangleVertexId >= 3)){
    stringstream msg;
    msg << "[Triangle] Triangle vertex Id #" << triangleVertexId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);

    return NULL;
  }

  if((mesh_)&&(vertexIds_))
    return mesh_->getVertex((*vertexIds_)[triangleVertexId]);
  else 
    return NULL;

  return 0;
}

int Triangle::getVertexId(const int &triangleVertexId, int &vertexId) const{

  if(!vertexIds_) return -1;

  if((triangleVertexId < 0)||(triangleVertexId >= 3)){
    stringstream msg;
    msg << "[Triangle] Triangle vertex Id #" << triangleVertexId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);

    return -2;
  }

  vertexId = (*vertexIds_)[triangleVertexId];

  return 0;
}

// TriangleMesh implementation
TriangleMesh::TriangleMesh(){

}

TriangleMesh::~TriangleMesh(){}

TriangleMesh::TriangleMesh(const vector<vector<double> > &pointList,
  const vector<vector<int> > &triangleList){

  set(pointList, triangleList);
}

TriangleMesh::TriangleMesh(const TriangleMesh &other){

  if(&other != this){

    set(*(other.getPointList()), *(other.getTriangleList()));
  }
}

TriangleMesh& TriangleMesh::operator=(const TriangleMesh &other){

  if(&other != this){
   
    set(*(other.getPointList()), *(other.getTriangleList()));
  }

  return *this;
}

void TriangleMesh::operator<<(istream &f){

  DebugTimer timer;

  int vertexNumber, faceNumber;
  string keyword;

  f >> keyword;
 
  if(keyword != "OFF")
    dMsg(cerr, 
      "[TriangleMesh] Input file does not seem to be a valid *off file.\n", 0);

  f >> vertexNumber;
  f >> faceNumber;
  f >> keyword;

  rawPointList_.resize(vertexNumber);
  rawTriangleList_.resize(faceNumber);

  vertices_.resize(vertexNumber);
  triangles_.resize(faceNumber); 
   
  // allocate memory for vertex locations
  for(int i = 0; i < rawPointList_.size(); i++){
    rawPointList_[i].resize(3);
    vertices_[i].point_ = &(rawPointList_[i]);
  }

  for(int i = 0; i < rawTriangleList_.size(); i++){
    rawTriangleList_[i].resize(3);
    triangles_[i].vertexIds_ = &(rawTriangleList_[i]);
  }

  for(int i = 0; i < vertexNumber; i++){
    vertices_[i].mesh_ = this;
    vertices_[i] << f;
  }

  for(int i = 0; i < faceNumber; i++){
    triangles_[i].mesh_ = this;
    triangles_[i] << f;
  }


  updateConnectivity();
  updateBoundingBox();

  stringstream msg;

  msg << "[TriangleMesh] Surface read in " 
    << timer.getElapsedTime()
    << " s. ("
    << vertices_.size() << " v, "
    << edges_.size() << " e, "
    << triangles_.size() << " t)."
    << endl; 

  dMsg(cout, msg.str(), 1);
}

void TriangleMesh::operator>>(ostream &f) const{
  
  f << "OFF" << endl;
  f << vertices_.size() << " " << triangles_.size() << " 0" << endl;

  for(int i = 0; i < vertices_.size(); i++)
    vertices_[i] >> f;

  for(int i = 0; i < triangles_.size(); i++)
    triangles_[i] >> f;

}

int TriangleMesh::removeIsolatedVertices(){

  vector<int>             nonIsolatedVertexIds(vertices_.size(), -1);
  vector<vector<double> > nonIsolatedPoints;

  const Vertex *v = NULL;
  const vector<double> *p = NULL;

  const Triangle *t = NULL;
  vector<int> triangleIds(3);
  
  for(int i = 0; i < triangles_.size(); i++){
    t = getTriangle(i);

    for(int j = 0; j < triangleIds.size(); j++){

      t->getVertexId(j, triangleIds[j]);
      if(nonIsolatedVertexIds[triangleIds[j]] == -1){

        nonIsolatedVertexIds[triangleIds[j]] = nonIsolatedPoints.size();
        v = getVertex(triangleIds[j]);
        p = v->getPoint();

        nonIsolatedPoints.push_back(*p);
      }
    }
  }

  vector<vector<int> > triangleList;

  for(int i = 0; i < triangles_.size(); i++){
    t = getTriangle(i);

    for(int j = 0; j < triangleIds.size(); j++){
      t->getVertexId(j, triangleIds[j]);
      triangleIds[j] = nonIsolatedVertexIds[triangleIds[j]];
    }

    triangleList.push_back(triangleIds);
  }

  set(nonIsolatedPoints, triangleList);

  return 0;
}

double TriangleMesh::rescaleBoundingBox(const double &sizeX, const double &sizeY,
  const double &sizeZ){

  double globalScale = 1;
  double maxEdge = bBoxX_.second - bBoxX_.first;

  globalScale = sizeX/maxEdge;

  if(bBoxY_.second - bBoxY_.first > maxEdge){
    maxEdge = bBoxY_.second - bBoxY_.first;
    globalScale = sizeY/maxEdge;
  }

  if(bBoxZ_.second - bBoxZ_.first > maxEdge){
    maxEdge = bBoxZ_.second - bBoxZ_.first;
    globalScale = sizeZ/maxEdge;
  }

  rescaleBoundingBox(globalScale);

  return globalScale;
}

int TriangleMesh::rescaleBoundingBox(const double &ratio){

  for(int i = 0; i < rawPointList_.size(); i++){
    for(int j = 0; j < rawPointList_[i].size(); j++){
      rawPointList_[i][j] *= ratio;
    }
  }

  updateBoundingBox();

  return 0;
}

int TriangleMesh::setBarycenter(const double &x, const double &y, 
  const double &z){

  if(!rawPointList_.size()) return -1;

  updateBoundingBox();

  vector<double> bary(3);

  bary[0] = (bBoxX_.second + bBoxX_.first)/2.0;
  bary[1] = (bBoxY_.second + bBoxY_.first)/2.0;
  bary[2] = (bBoxZ_.second + bBoxZ_.first)/2.0;

  for(int i = 0; i < rawPointList_.size(); i++){
    rawPointList_[i][0] += x - bary[0];
    rawPointList_[i][1] += y - bary[1];
    rawPointList_[i][2] += z - bary[2];
  }

  return 0;  
}

int TriangleMesh::updateBoundingBox(){

  if(!rawPointList_.size()) return -1;

  bBoxX_.first = bBoxX_.second = bBoxY_.first = bBoxY_.second 
    = bBoxZ_.first = bBoxZ_.second = -1;

  // base the loop on the faces to avoid isolated vertices if any

  vector<bool>  visitedVertices(rawPointList_.size(), false);

  int pId = 0;

  for(int i = 0; i < rawTriangleList_.size(); i++){

    for(int vId = 0; vId < 3; vId++){
      
      pId = rawTriangleList_[i][vId];

      if(!visitedVertices[pId]){

        visitedVertices[pId] = true;
    
        if(((!i)&&(!vId))||(rawPointList_[pId][0] < bBoxX_.first))
          bBoxX_.first = rawPointList_[pId][0];
        if(((!i)&&(!vId))||(rawPointList_[pId][0] > bBoxX_.second))
          bBoxX_.second = rawPointList_[pId][0];
    
        if(((!i)&&(!vId))||(rawPointList_[pId][1] < bBoxY_.first))
          bBoxY_.first = rawPointList_[pId][1];
        if(((!i)&&(vId))||(rawPointList_[pId][1] > bBoxY_.second))
          bBoxY_.second = rawPointList_[pId][1];
    
        if(((!i)&&(!vId))||(rawPointList_[pId][2] < bBoxZ_.first))
          bBoxZ_.first = rawPointList_[pId][2];
        if(((!i)&&(vId))||(rawPointList_[pId][2] > bBoxZ_.second))
          bBoxZ_.second = rawPointList_[pId][2];
      }
    }
  }

  return 0;
}

int TriangleMesh::updateConnectivity(){
  
  stringstream                  msg;

  map<std::pair<int, int>, int> edgeMap;

  for(int i = 0; i < triangles_.size(); i++){

    for(int j = 0; j < 3; j++){

      // add the triangle to its vertex
      vertices_[(*triangles_[i].vertexIds_)[j]].triangleIds_.push_back(i);

      pair<int, int> eInfo;
      if(!triangles_[i].vertexIds_) break;

      eInfo.first = (*triangles_[i].vertexIds_)[j];
      eInfo.second = (*triangles_[i].vertexIds_)[(j + 1)%3];

      map<pair<int, int>, int>::iterator eIt = edgeMap.find(eInfo);
      if(eIt == edgeMap.end()){
        eInfo.first = (*triangles_[i].vertexIds_)[(j + 1)%3];
        eInfo.second = (*triangles_[i].vertexIds_)[j];
        eIt = edgeMap.find(eInfo);
      }
      if(eIt == edgeMap.end()){
        // definitely not in there

        edgeMap[eInfo] = edges_.size();
          
        Edge e;
        e.vertexIds_[0] = (*triangles_[i].vertexIds_)[j];
        e.vertexIds_[1] = (*triangles_[i].vertexIds_)[(j + 1)%3];
        e.triangleIds_[0] = i;

        // update vertices
        vertices_[e.vertexIds_[0]].starVertexIds_.push_back(e.vertexIds_[1]);
        vertices_[e.vertexIds_[1]].starVertexIds_.push_back(e.vertexIds_[0]);
        vertices_[e.vertexIds_[0]].edgeIds_.push_back(edges_.size());
        vertices_[e.vertexIds_[1]].edgeIds_.push_back(edges_.size());

        // update triangles
        triangles_[i].edgeIds_[j] = edges_.size();
        
        e.mesh_ = this;
        edges_.push_back(e);
      }
      else{
        // the edge is already in there.

        // update triangles
        triangles_[i].edgeIds_[j] = eIt->second;
        
        //update the edge
        if(edges_[eIt->second].triangleIds_[1] != -1){
          msg << "[TriangleMesh] Warning, non-manifold edge #" << eIt->second
            << "flushes connection to triangle #" << 
              edges_[eIt->second].triangleIds_[1] << endl;
          dMsg(cerr, msg.str(), 10);
        }
        edges_[eIt->second].triangleIds_[1] = i;
      }
    }
  }

  // update that for convenience
  for(int i = 0; i < triangles_.size(); i++){
    for(int j = 0; j < 3; j++){
      int n = edges_[triangles_[i].edgeIds_[j]].triangleIds_[0];
      if(n == i){
        n = edges_[triangles_[i].edgeIds_[j]].triangleIds_[1];
      }

      triangles_[i].triangleIds_[j] = n;
    }
  }
 
  // update opposite edges
  int edgeId;
  for(int i = 0; i < vertices_.size(); i++){
    for(int j = 0; j < vertices_[i].triangleIds_.size(); j++){
      for(int k = 0; k < 3; k++){ 
        edgeId = triangles_[vertices_[i].triangleIds_[j]].edgeIds_[k];
        if((edges_[edgeId].vertexIds_[0] != i)
          &&(edges_[edgeId].vertexIds_[1] != i)){
          vertices_[i].linkEdgeIds_.push_back(edgeId);
        }
      }
    }
  }

  // now re-order vertices' star triangles and link edges for fast traversal.
  int           pivotVertexId = -1, edgePivot, visitedEdgeNb;
  vector<bool>  visitedEdges;
  vector<int>   orderedStarTriangles, orderedLinkEdges;
  Edge          *e = NULL, *nextEdge = NULL;
  
  for(int i = 0; i < vertices_.size(); i++){
    
    if(!vertices_[i].linkEdgeIds_.size()) break;

    orderedStarTriangles.resize(vertices_[i].triangleIds_.size());
    orderedLinkEdges.resize(vertices_[i].linkEdgeIds_.size());
    visitedEdges.resize(vertices_[i].linkEdgeIds_.size());
    for(int j = 0; j < visitedEdges.size(); j++) visitedEdges[j] = false;
    visitedEdgeNb = 0;
    pivotVertexId = -1;

    // start with a boundary vertex if any
    for(int j = 0; j < vertices_[i].linkEdgeIds_.size(); j++){
      if(vertices_[
        edges_[vertices_[i].linkEdgeIds_[j]].vertexIds_[0]].isOnBoundary()){
        
        pivotVertexId = edges_[vertices_[i].linkEdgeIds_[j]].vertexIds_[1];
        edgePivot = 1;
        e = &(edges_[vertices_[i].linkEdgeIds_[j]]);
        visitedEdges[j] = true;
        orderedLinkEdges[0] = vertices_[i].linkEdgeIds_[j];
        visitedEdgeNb++;
        break;
      }
      if(vertices_[
        edges_[vertices_[i].linkEdgeIds_[j]].vertexIds_[1]].isOnBoundary()){
        
        pivotVertexId = edges_[vertices_[i].linkEdgeIds_[j]].vertexIds_[0];
        edgePivot = 0;
        e = &(edges_[vertices_[i].linkEdgeIds_[j]]);
        visitedEdges[j] = true;
        orderedLinkEdges[0] = vertices_[i].linkEdgeIds_[j];
        visitedEdgeNb++;
        break;
      }
    }

    if(pivotVertexId == -1){
      // no boundary vertex
      e = &(edges_[vertices_[i].linkEdgeIds_[0]]);
      visitedEdges[0] = true;
      orderedLinkEdges[0] = vertices_[i].linkEdgeIds_[0];
      edgePivot = 0;
      visitedEdgeNb++;
    }

    while(visitedEdgeNb < orderedLinkEdges.size()){
     
      e->getVertexId(edgePivot, pivotVertexId);

      // find the next edge to add to the list
      e = NULL;
      for(int j = 0; j < vertices_[i].linkEdgeIds_.size(); j++){
        if(!visitedEdges[j]){
          nextEdge = &(edges_[vertices_[i].linkEdgeIds_[j]]);
          if((nextEdge->vertexIds_[0] == pivotVertexId)
            ||(nextEdge->vertexIds_[1] == pivotVertexId)){
          
            orderedLinkEdges[visitedEdgeNb] = vertices_[i].linkEdgeIds_[j];
            visitedEdges[j] = true;
            visitedEdgeNb++;
          
            e = nextEdge;

            if(nextEdge->vertexIds_[0] == pivotVertexId) edgePivot = 1;
            if(nextEdge->vertexIds_[1] == pivotVertexId) edgePivot = 0;
          }
        }
      }
      if((!e)&&(visitedEdgeNb < orderedLinkEdges.size())){
        // we didn't find an edge common to the pivot vertex and still we didn't
        // visit everything (we reached the boundary of the surface)

        // let's try to catch one that touches the boundary
        // (there has to be one)
        for(int j = 0; j < vertices_[i].linkEdgeIds_.size(); j++){
          if(!visitedEdges[j]){
            if(vertices_[
              edges_[vertices_[i].linkEdgeIds_[
                j]].vertexIds_[0]].isOnBoundary()){
        
              pivotVertexId = edges_[
                vertices_[i].linkEdgeIds_[j]].vertexIds_[1];
              edgePivot = 1;
              e = &(edges_[vertices_[i].linkEdgeIds_[j]]);
              visitedEdges[j] = true;
              orderedLinkEdges[visitedEdgeNb] = vertices_[i].linkEdgeIds_[j];
              visitedEdgeNb++;
              break;
            }
            if(vertices_[
              edges_[vertices_[i].linkEdgeIds_[
                j]].vertexIds_[1]].isOnBoundary()){
        
              pivotVertexId = edges_[
                vertices_[i].linkEdgeIds_[j]].vertexIds_[0];
              edgePivot = 0;
              e = &(edges_[vertices_[i].linkEdgeIds_[j]]);
              visitedEdges[j] = true;
              orderedLinkEdges[visitedEdgeNb] = vertices_[i].linkEdgeIds_[j];
              visitedEdgeNb++;
              break;
            }
          }
        }
        if(!e){
          stringstream msg;
          msg 
            << "[TriangleMesh] UpdateConnectivity: no next bEdge on vertex #" 
            << i << "!"  << endl;
          dMsg(cerr, msg.str(), 1);
          return 0;
        }
      }
    }

    vertices_[i].linkEdgeIds_ = orderedLinkEdges;

    visitedEdgeNb = 0;
    // re-order star triangles accordingly
    for(int j = 0; j < vertices_[i].linkEdgeIds_.size(); j++){
      for(int k = 0; k < vertices_[i].triangleIds_.size(); k++){
        if(triangles_[vertices_[i].triangleIds_[k]].edgeIds_[0] == 
          vertices_[i].linkEdgeIds_[j]){

          orderedStarTriangles[visitedEdgeNb] = 
            vertices_[i].triangleIds_[k];
          visitedEdgeNb++;
          break;
        }
        else if(triangles_[vertices_[i].triangleIds_[k]].edgeIds_[1] == 
          vertices_[i].linkEdgeIds_[j]){

          orderedStarTriangles[visitedEdgeNb] = 
            vertices_[i].triangleIds_[k];
          visitedEdgeNb++;
          break;
        }
        else if(triangles_[vertices_[i].triangleIds_[k]].edgeIds_[2] == 
          vertices_[i].linkEdgeIds_[j]){

          orderedStarTriangles[visitedEdgeNb] = 
            vertices_[i].triangleIds_[k];
          visitedEdgeNb++;
          break;
        }
      }
    }
    vertices_[i].triangleIds_ = orderedStarTriangles;
  }

  // now re-order vertices' star edges and star vertices
  int           pivotTriangleId;

  vector<int>   orderedEdges;
  vector<int>   orderedVertices;
  for(int i = 0; i < vertices_.size(); i++){
  
    if(!vertices_[i].edgeIds_.size()) break;

    orderedEdges.resize(vertices_[i].edgeIds_.size());
    orderedVertices.resize(vertices_[i].starVertexIds_.size());
    visitedEdges.resize(vertices_[i].edgeIds_.size());
    for(int j = 0; j < visitedEdges.size(); j++) visitedEdges[j] = false;

    pivotTriangleId = -1;
    visitedEdgeNb = 0;

    // start with a boundary edge if any
    for(int j = 0; j < vertices_[i].edgeIds_.size(); j++){
      if(edges_[vertices_[i].edgeIds_[j]].isOnBoundary()){
        e =  &(edges_[vertices_[i].edgeIds_[j]]);
        orderedEdges[visitedEdgeNb] = vertices_[i].edgeIds_[j];
        visitedEdgeNb++;
        if(e->triangleIds_[0] == -1){
          pivotTriangleId = e->triangleIds_[1];
        }
        else if(e->triangleIds_[1] == -1){
          pivotTriangleId = e->triangleIds_[0];
        }
        visitedEdges[j] =true;
        break;
      }
    }

    if(pivotTriangleId == -1){
      // non boundary vertex
      e = &(edges_[vertices_[i].edgeIds_[0]]);
      orderedEdges[visitedEdgeNb] = vertices_[i].edgeIds_[0];
      visitedEdgeNb++;
      pivotTriangleId = e->triangleIds_[0];
      visitedEdges[0] = true;
    }

    // as is this code should handle boundary (and re-loop starting at a
    // boundary)
    bool complete = false;
    while(visitedEdgeNb < orderedEdges.size()){
      complete = false;
      for(int j = 0; j < vertices_[i].edgeIds_.size(); j++){
        if(!visitedEdges[j]){
          // it's not the last visited edge
          e = &(edges_[vertices_[i].edgeIds_[j]]);
          if(e->triangleIds_[0] == pivotTriangleId){
            orderedEdges[visitedEdgeNb] = vertices_[i].edgeIds_[j];
            visitedEdgeNb++;
            pivotTriangleId = e->triangleIds_[1];
            visitedEdges[j] = true;
            complete = true;
            break;
          }
          if(e->triangleIds_[1] == pivotTriangleId){
            orderedEdges[visitedEdgeNb] = vertices_[i].edgeIds_[j];
            visitedEdgeNb++;
            pivotTriangleId = e->triangleIds_[0];
            visitedEdges[j] = true;
            complete = true;
            break;
          }
        }
      }
      if(!complete){
        // the mesh is not manifold
        // we looped over the vertex's edges and did not find any new candidate
        // the while loop will continue forever if we don't stop.
        break;
      }
    }
    vertices_[i].edgeIds_ = orderedEdges;

    // now make the vertices
    visitedEdgeNb = 0;

    for(int j = 0; j < vertices_[i].edgeIds_.size(); j++){
      if(edges_[vertices_[i].edgeIds_[j]].vertexIds_[0] == i)
        vertices_[i].starVertexIds_[visitedEdgeNb] 
          = edges_[vertices_[i].edgeIds_[j]].vertexIds_[1];
      else if(edges_[vertices_[i].edgeIds_[j]].vertexIds_[1] == i)
        vertices_[i].starVertexIds_[visitedEdgeNb] 
          = edges_[vertices_[i].edgeIds_[j]].vertexIds_[0];
      visitedEdgeNb++;
    }
  }


  return 0;
}

int TriangleMesh::exportToVRML(ostream &f, 
  const vector<pair<double, double> > &uvList, 
  const string &textureFileString) const{
  
  f << "\t\tShape{" << endl;

  // default appearance
  f << "\t\t\tappearance Appearance{" << endl;
  f << "\t\t\t\tmaterial Material{" << endl;
  f << "\t\t\t\t\tdiffuseColor 0.7 0.7 0.7" << endl;
  f << "\t\t\t\t\tspecularColor 0.7 0.7 0.7" << endl;
  f << "\t\t\t\t\tshininess 0.9" << endl;
  f << "\t\t\t\t}" << endl;

  if(uvList.size() == vertices_.size()){
    
    f << "\t\t\t\ttexture ImageTexture{" << endl;
    f << "\t\t\t\t\turl \"" << textureFileString << "\"" << endl;
    f << "\t\t\t\t\trepeatS TRUE" << endl;
    f << "\t\t\t\t\trepeatT TRUE" << endl;
    f << "\t\t\t\t}" << endl;
    
  }
  f << "\t\t\t}" << endl;

  // geometry
  f << "\t\t\tgeometry IndexedFaceSet{" << endl;
  f << "\t\t\t\tcoord Coordinate{" << endl;
  f << "\t\t\t\t\tpoint[" << endl;
  for(int i = 0; i < vertices_.size(); i++)
    vertices_[i].exportToVRML(f);
  f << "\t\t\t\t\t]" << endl;
  f << "\t\t\t\t}" << endl;

  f << "\t\t\t\tcoordIndex[" << endl;
  for(int i = 0; i < triangles_.size(); i++)
    triangles_[i].exportToVRML(f);
  f << "\t\t\t\t]" << endl;

  if(uvList.size() == vertices_.size()){

    f << "\t\t\t\ttexCoord TextureCoordinate{" << endl;
    f << "\t\t\t\t\tpoint[" << endl;
    for(int i = 0; i < uvList.size(); i++)
      f << "\t\t\t\t\t\t " << uvList[i].first 
        << " " << uvList[i].second << endl;
    f << "\t\t\t\t\t]" << endl;
    f << "\t\t\t\t}" << endl;

    f << "\t\t\t\ttexCoordIndex[" << endl;
    for(int i = 0; i < triangles_.size(); i++)
      triangles_[i].exportToVRML(f);
    f << "\t\t\t\t]" << endl;
  }

  f << "\t\t\t}" << endl;
  
  f << "\t\t}" << endl;

  return 0;
}

const Edge* TriangleMesh::getEdge(const int &edgeId) const{

  if((edgeId < 0)||(edgeId >= edges_.size())){
    stringstream msg;
    msg << "[TriangleMesh] Edge Id #" << edgeId << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return NULL;
  }

  return (const Edge *)(&edges_[edgeId]);
}

double TriangleMesh::getHausdorffApproximationTo(
  const TriangleMesh &other) const{

  DebugTimer timer;

  int             measureNb = 0;

  double          currentDistance, minVertexDistance = 1, maxDistance = 1,
                  average = 0,
                  d0, d1, d2, triangleDistance;
  vector<double>  p(3), 
                  p0(3), p1(3), p2(3),
                  projectedPoint(3);

  vector<bool>    visitedVertices(rawPointList_.size(), false);

  // base the loop on the faces to avoid isolated vertices
  for(int i = 0; i < rawTriangleList_.size(); i++){
  
    for(int vId = 0; vId < 3; vId++){

      if(!visitedVertices[rawTriangleList_[i][vId]]){

        visitedVertices[rawTriangleList_[i][vId]] = true;
        p = rawPointList_[rawTriangleList_[i][vId]];

        for(int j = 0; j < other.triangles_.size(); j++){

          currentDistance = 0;

          // get the points
          p0 = other.rawPointList_[other.rawTriangleList_[j][0]];
          p1 = other.rawPointList_[other.rawTriangleList_[j][1]];
          p2 = other.rawPointList_[other.rawTriangleList_[j][2]];
     
          d0 = 0;
          for(int k = 0; k < 3; k++) d0 += (p[k] - p0[k])*(p[k] - p0[k]);
          d0 = sqrt(d0);
          currentDistance = d0;

          d1 = 0;
          for(int k = 0; k < 3; k++) d1 += (p[k] - p1[k])*(p[k] - p1[k]);
          d1 = sqrt(d1);

          if(d1 < currentDistance) currentDistance = d1;

          d2 = 0;
          for(int k = 0; k < 3; k++) d2 += (p[k] - p2[k])*(p[k] - p2[k]);
          d2 = sqrt(d2);
          if(d2 < currentDistance) currentDistance = d2;

          if(((!j)
            ||(currentDistance < minVertexDistance))
            &&(!isnan(currentDistance)))
            minVertexDistance = currentDistance;
        }
      }
    
   
      average += minVertexDistance;
      measureNb++;

      if((((!i)&&(!vId))
        ||(minVertexDistance > maxDistance))&&(!isnan(minVertexDistance)))
        maxDistance = minVertexDistance;
    }
  }

  average /= (double (measureNb));

  stringstream msg;

  msg << "[TriangleMesh] Hausdorff approximation: " << maxDistance 
    << " (computed in  " 
    << timer.getElapsedTime()
    << " s.)" << endl;
  msg << "[TriangleMesh] Average: " << average << "" << endl;


  dMsg(cout, msg.str(), 1);

  return maxDistance;
}


double TriangleMesh::getHausdorffTo(const TriangleMesh &other) const{

  DebugTimer timer;

  int             measureNb = 0;

  double          currentDistance, minVertexDistance = 1, maxDistance = 1,
                  average = 0,
                  d0, d1, d2, triangleDistance;
  vector<double>  p(3), 
                  p0(3), p1(3), p2(3),
                  projectedPoint(3);

  vector<bool>    visitedVertices(rawPointList_.size(), false);
  int pId = 0;

  // base the loop on the faces to avoid potential isolated vertices
  for(int i = 0; i < rawTriangleList_.size(); i++){
  
    for(int vId = 0; vId < 3; vId++){

      pId = rawTriangleList_[i][vId];

      if(!visitedVertices[pId]){

        visitedVertices[pId] = true;

        p = rawPointList_[pId];

        for(int j = 0; j < other.triangles_.size(); j++){

          currentDistance = 0;

          // get the points
          p0 = other.rawPointList_[other.rawTriangleList_[j][0]];
          p1 = other.rawPointList_[other.rawTriangleList_[j][1]];
          p2 = other.rawPointList_[other.rawTriangleList_[j][2]];
     
          d0 = 0;
          for(int k = 0; k < 3; k++) d0 += (p[k] - p0[k])*(p[k] - p0[k]);
          d0 = sqrt(d0);
          currentDistance = d0;

          d1 = 0;
          for(int k = 0; k < 3; k++) d1 += (p[k] - p1[k])*(p[k] - p1[k]);
          d1 = sqrt(d1);

          if(d1 < currentDistance) currentDistance = d1;

          d2 = 0;
          for(int k = 0; k < 3; k++) d2 += (p[k] - p2[k])*(p[k] - p2[k]);
          d2 = sqrt(d2);
          if(d2 < currentDistance) currentDistance = d2;

          other.projectPointOnTrianglePlane(j, &p, &projectedPoint);

          if(other.isPointInTriangle(j, &projectedPoint)){

            triangleDistance = 0;

            for(int k = 0; k < 3; k++) 
              triangleDistance += 
                (projectedPoint[k] - p[k])*(projectedPoint[k] - p[k]);
            triangleDistance = sqrt(triangleDistance);

            if(triangleDistance > currentDistance){
              stringstream msg;
              msg << 
          "[TriangleMesh] getHausdorff: projected point further than vertices!"
              << endl;
              dMsg(cerr, msg.str(), 10);
            }
            else{
              currentDistance = triangleDistance;
            }        
          }

          if(((!j)
            ||(currentDistance < minVertexDistance))
            &&(!isnan(currentDistance))){
            minVertexDistance = currentDistance;
          }
        }
      }
    
   
      average += minVertexDistance;
      measureNb++;

      if((((!i)&&(!vId))||(minVertexDistance > maxDistance))
        &&(!isnan(minVertexDistance)))
        maxDistance = minVertexDistance;
    }
  }

  average /= ((double) measureNb);

  stringstream msg;

  msg << "[TriangleMesh] Hausdorff distance: " << maxDistance 
    << " (computed in  " 
    << timer.getElapsedTime()
    << " s.)" << endl;
  msg << "[TriangleMesh] Average: " << average << "" << endl;

  dMsg(cout, msg.str(), 1);

  return maxDistance;
}

const Triangle* TriangleMesh::getTriangle(const int &triangleId) const{

  if((triangleId < 0)||(triangleId >= triangles_.size())){
    stringstream msg;
    msg << "[TriangleMesh] Triangle Id #" << triangleId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return NULL;
  }

  return (const Triangle *) &triangles_[triangleId];
}

int TriangleMesh::getVertexId(const int &triangleId, const int &triangleEdgeId, 
  const int &edgeVertexId, int &vertexId) const{

  if((triangleId < 0)||(triangleId >= triangles_.size())){
    stringstream msg;
    msg << "[TriangleMesh] Triangle Id #" << triangleId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return -1;
  }
  if((triangleEdgeId < 0)||(triangleEdgeId > 2)){
    stringstream msg;
    msg << "[TriangleMesh] Edge Id #" << triangleEdgeId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return -2;
  }

  const Edge *e = triangles_[triangleId].getEdge(triangleEdgeId);
  if(e) e->getVertexId(edgeVertexId, vertexId);

  return 0;
}

int TriangleMesh::set(const vector<vector<double> > &pointList, 
  const vector<vector<int> > &triangleList){

  DebugTimer timer;

  rawPointList_.resize(pointList.size());
  rawTriangleList_.resize(triangleList.size());

  vertices_.resize(pointList.size());
  edges_.clear();
  triangles_.resize(triangleList.size());

  // allocate memory for vertex locations
  for(int i = 0; i < pointList.size(); i++){
    rawPointList_[i].resize(3);
    rawPointList_[i] = pointList[i];
    vertices_[i].point_ = &(rawPointList_[i]);
    vertices_[i].mesh_ = this;
  }

  for(int i = 0; i < triangleList.size(); i++){
    rawTriangleList_[i].resize(3);
    rawTriangleList_[i] = triangleList[i];
    triangles_[i].vertexIds_ = &(rawTriangleList_[i]);
    triangles_[i].mesh_ = this;
  }

  updateConnectivity();
  updateBoundingBox();

  stringstream msg;

  msg << "[TriangleMesh] Surface set in " 
    << timer.getElapsedTime()
    << " s. ("
    << vertices_.size() << " v, "
    << edges_.size() << " e, "
    << triangles_.size() << " t)."
    << endl; 

  dMsg(cout, msg.str(), 1);

  return 0;  
}
