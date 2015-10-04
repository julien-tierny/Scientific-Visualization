/*
 * file:                TriangleMesh.h
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

#ifndef                 _TRIANGLE_MESH_H
#define                 _TRIANGLE_MESH_H

#include                <Debug.h>

#include                <iostream>
#include                <map>
#include                <sstream>
#include                <string>
#include                <vector>

#include                <math.h>

using namespace std;

class Edge;

class Triangle;

class TriangleMesh;

class Vertex : public Debug{

  public:
    
    friend class TriangleMesh;

    Vertex();

    Vertex(const vector<double> &point);

    Vertex(const Vertex &other);

    Vertex(const Vertex &other, TriangleMesh *mesh);

    ~Vertex();

    Vertex& operator=(const Vertex &other);

    void operator<<(istream &f);

    void operator>>(ostream &f) const;

    int exportToVRML(ostream &f) const;

    const Edge* getEdge(const int &vertexEdgeId) const;

    int getEdgeId(const int &vertexEdgeId, int &edgeId) const;

    const Edge* getLinkEdge(const int &linkEdgeId) const;

    int getLinkEdgeId(const int &linkEdgeId, int &edgeId) const;

    inline int getNumberOfEdges() const { return edgeIds_.size();};

    inline int getNumberOfBoundaryEdges() const;

    inline int getNumberOfLinkEdges() const { return linkEdgeIds_.size();};

    inline int getNumberOfStarVertices() const { return starVertexIds_.size();};

    inline int getNumberOfTriangles() const { return triangleIds_.size();};

    inline const vector<double>* getPoint() const { return point_;};

    inline int getPoint(vector<double> &p) const {
      if(!point_) return -1;
      p.resize(point_->size());
      for(int i = 0; i < point_->size(); i++) p[i] = (*point_)[i];
      return 0;};

    inline const Vertex* getStarVertex(const int &starVertexId) const;

    inline int getStarVertexId(const int &starVertexId, int &vertexId) const;

    inline const Triangle* getTriangle(const int &vertexTriangleId) const;

    int getTriangleId(const int &vertexTriangleId, int &triangleId) const;

    inline bool isOnBoundary() const {
      return (getNumberOfBoundaryEdges() != 0);};


  protected:

    vector<double>                *point_;
    vector<int>                   starVertexIds_;
    vector<int>                   edgeIds_, linkEdgeIds_;
    vector<int>                   triangleIds_;
    TriangleMesh                  *mesh_;      
};

class Edge : public Debug{
  
  public:

    friend class TriangleMesh;
    
    Edge();

    Edge(const Edge &other);

    Edge(const Edge &other, TriangleMesh *mesh);

    ~Edge();

    Edge& operator=(const Edge &other);

    void operator<<(istream &f);

    void operator>>(ostream &f) const;

    const Triangle* getTriangle(const int &edgeTriangleId) const;

    int getTriangleId(const int &edgeTriangleId, int &triangleId) const;

    const Vertex* getVertex(const int &edgeVertexId) const;

    int getVertexId(const int &edgeVertexId, int &vertexId) const;

    inline bool isOnBoundary() const{ 
      return ((triangleIds_[0] == -1)||(triangleIds_[1] == -1));};

  protected:
    
    int                           vertexIds_[2];
    int                           triangleIds_[2];
    TriangleMesh                  *mesh_;
};

class Triangle : public Debug{

  public:
    
    friend class TriangleMesh;

    Triangle();

    Triangle(const vector<int> &vertexList);

    Triangle(const Triangle &other);

    Triangle(const Triangle &other, TriangleMesh *mesh);

    ~Triangle();

    Triangle& operator=(const Triangle &other);

    void operator<<(istream &f);

    void operator>>(ostream &f) const;

    int exportToVRML(ostream &f) const;

    const Edge* getEdge(const int &triangleEdgeId) const;

    int getEdgeId(const int &triangleEdgeId, int &edgeId) const;

    inline int getPivotEdgeId(const int &triangleId, int &pivotEdgeId) const;

    const Triangle* getTriangle(const int &triangleTriangleId) const;

    int getTriangleId(const int &triangeTriangleId, int &triangleId) const;

    const Vertex* getVertex(const int &triangleVertexId) const;

    int getVertexId(const int &triangleVertexId, int &vertexId) const;

    inline const vector<int>* getVertexIds() const { return vertexIds_;};

  protected:

    vector<int>                   *vertexIds_;
    int                           edgeIds_[3];
    int                           triangleIds_[3];
    TriangleMesh                  *mesh_;
};

class TriangleMesh : public Debug {

  public:

    // standard accessors
    TriangleMesh();

    TriangleMesh(const vector<vector<double> > &pointList, 
      const vector<vector<int> > &triangleList);

    TriangleMesh(const TriangleMesh &other);

    ~TriangleMesh();

    TriangleMesh& operator=(const TriangleMesh &other);

    void operator<<(istream &f);

    void operator>>(ostream &f) const;

    int exportToVRML(ostream &f, const vector<pair<double, double> > &uvList,
      const string &textureFileString) const;

    inline int getBoundingBox(pair<double, double> &bBoxX, 
      pair<double, double> &bBoxY, pair<double, double> &bBoxZ) const{
      bBoxX = bBoxX_; bBoxY = bBoxY_; bBoxZ = bBoxZ_;
      return 0;  
    }
 
    const Edge* getEdge(const int &edgeId) const;

    double getHausdorffApproximationTo(const TriangleMesh &other) const;

    double getHausdorffTo(const TriangleMesh &other) const;

    inline int getNumberOfEdges() const {return edges_.size();};

    inline int getNumberOfTriangles() const{ return triangles_.size();};

    inline int getNumberOfVertices() const{ return vertices_.size();};

    inline int getPoint(const int &triangleId, const int &triangleEdgeId,
      const double &edgeTime, vector<double> *p) const;

    inline int getPoint(const int &vertexId, vector<double> *p) const{
      if((vertexId < 0)||(vertexId >= vertices_.size())) return -1;
      if(!p) return -2;
      p->resize(3);
      for(int i = 0; i < 3; i++) (*p)[i] = rawPointList_[vertexId][i];
      return 0;
    }

    inline const vector<vector<double> >* getPointList() const{
      return (const vector<vector<double> >*)
        &rawPointList_;};

    const Triangle* getTriangle(const int &triangleId) const;

    inline int getTriangleArea(const int &triangleId,
      double &area) const;

    inline const vector<vector<int> >* getTriangleList() const{
      return (const vector<vector<int> >*)
        &rawTriangleList_;};

    inline int getTriangleNormal(const int &triangleId, 
      vector<double> *normal) const;

    inline const Vertex* getVertex(const int &vertexId) const;

    int getVertexId(const int &triangleId, const int &triangleEdgeId, 
      const int &edgeVertexId, int &vertexId) const;

    inline bool isPointInTriangle(const int &triangleId, 
      const vector<double> *p) const;

    inline int projectPointOnTrianglePlane(const int &triangleId, 
      const vector<double> *point, vector<double> *projectedPoint) const;

    int removeIsolatedVertices();

    double rescaleBoundingBox(const double &sizeX, const double &sizeY,
      const double &sizeZ);

    int rescaleBoundingBox(const double &ratio);

    int set(const vector<vector<double> > &pointList, 
      const vector<vector<int> > &triangleList);

    int setBarycenter(const double &x, const double &y, const double &z);

    inline int setVertexPoint(const int &vertexId, const vector<double> &p);

  protected:

    int updateBoundingBox();
    
    int updateConnectivity();

    pair<double, double>          bBoxX_, bBoxY_, bBoxZ_;
    vector<vector<double> >       rawPointList_;
    vector<vector<int> >          rawTriangleList_;
    vector<Vertex>                vertices_;
    vector<Edge>                  edges_;
    vector<Triangle>              triangles_;
};

inline int Vertex::getNumberOfBoundaryEdges() const{

  int         boundaryEdgeNb = 0;
  const Edge  *e;

  if(!mesh_) return -1;

  for(int i = 0; i < edgeIds_.size(); i++){
    e = mesh_->getEdge(edgeIds_[i]);
    if(e){
      if(e->isOnBoundary()) boundaryEdgeNb++;
    }
  }

  return boundaryEdgeNb;
}

inline const Vertex* Vertex::getStarVertex(const int &starVertexId) const{

  if((starVertexId < 0)||(starVertexId >= starVertexIds_.size())){
    stringstream msg;
    msg << "[Vertex] Star vertex Id #" << starVertexId << " out of range."
      << endl;
    dMsg(cerr, msg.str(), 10);
    return NULL;
  }

  if(mesh_) return mesh_->getVertex(starVertexIds_[starVertexId]);
  else return NULL;
}

inline int Vertex::getStarVertexId(const int &starVertexId, 
  int &vertexId) const{

  if((starVertexId < 0)||(starVertexId >= starVertexIds_.size())){
    stringstream msg;
    msg << "[Vertex] Star vertex Id #" << starVertexId << " out of range."
      << endl;
    dMsg(cerr, msg.str(), 10);
    return -1;
  }

  vertexId = starVertexIds_[starVertexId];

  return 0;
}

inline const Triangle* Vertex::getTriangle(const int &vertexTriangleId) const{
  
  if((vertexTriangleId < 0)||(vertexTriangleId >= triangleIds_.size())){
    stringstream msg;
    msg << "[Vertex] Vertex triangle Id #" << vertexTriangleId 
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return NULL;
  }

  if(mesh_) return mesh_->getTriangle(triangleIds_[vertexTriangleId]);
  else return NULL;
}


inline int Triangle::getPivotEdgeId(const int &triangleId, 
  int &pivotEdgeId) const{

  for(int i = 0; i < 3; i++){
    if(triangleIds_[i] == triangleId){
      pivotEdgeId = i;
      return 0;
    }
  }

  return -1;
}

inline int TriangleMesh::getPoint(const int &triangleId, 
  const int &triangleEdgeId, const double &edgeTime, vector<double> *p) const{

  if((triangleId < 0)||(triangleId >= triangles_.size())) return -1;

  if((triangleEdgeId < 0)||(triangleEdgeId >= 3)) return -2;

  if((edgeTime < 0)||(edgeTime > 1)) return -3;

  if((!p)||(p->size() < 3)) return -4;

  for(int i = 0; i < 3; i++)
    (*p)[i] = (1 - edgeTime)*rawPointList_[
      edges_[triangles_[triangleId].edgeIds_[triangleEdgeId]].vertexIds_[0]][i]
      + (edgeTime)*rawPointList_[
      edges_[triangles_[triangleId].edgeIds_[triangleEdgeId]].vertexIds_[1]][i];

  return 0;
}

inline int TriangleMesh::getTriangleArea(const int &triangleId, 
  double &area) const{
  
  if((triangleId < 0)||(triangleId >= triangles_.size())) return -1;

  vector<double> normal(3);

  getTriangleNormal(triangleId, &normal);

  double nLength = 0;
  for(int i = 0; i < 3; i++) nLength += normal[i]*normal[i];
  nLength = sqrt(nLength);

  area = nLength/2;

  return 0;
 
}

inline int TriangleMesh::getTriangleNormal(const int &triangleId, 
  vector<double> *normal) const{

  if((triangleId < 0)||(triangleId >= triangles_.size())) return -1;

  if(!normal) return -2;

  normal->resize(3);

  vector<double>  e0(3), e1(3);

  for(int i = 0; i < 3; i++)
    e0[i] = rawPointList_[rawTriangleList_[triangleId][1]][i]
      - rawPointList_[rawTriangleList_[triangleId][0]][i];

  for(int i = 0; i < 3; i++)
    e1[i] = rawPointList_[rawTriangleList_[triangleId][2]][i]
      - rawPointList_[rawTriangleList_[triangleId][0]][i];

  // cross product
  for(int i = 0; i < 3; i++)
    (*normal)[i] = e0[(i+1)%3]*e1[(i+2)%3] - e0[(i+2)%3]*e1[(i+1)%3];
  
  return 0;
}

inline const Vertex* TriangleMesh::getVertex(const int &vertexId) const{

  if((vertexId < 0)||(vertexId >= vertices_.size())){
    stringstream msg;
    msg << "[TriangleMesh] Vertex Id #" << vertexId << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return NULL;
  }

  return (const Vertex *) &vertices_[vertexId];
}



inline bool TriangleMesh::isPointInTriangle(const int &triangleId, 
  const vector<double> *p) const{
  
  if((triangleId < 0)||(triangleId >= triangles_.size())) return false;

  if((!p)||(p->size() != 3)) return false;
 
  vector<double>    n(3);
  double            nLength = 0;
  getTriangleNormal(triangleId, &n);
  for(int i = 0; i < 3; i++) nLength += n[i]*n[i];
  nLength = sqrt(nLength);

  vector<double>    e0(3), e1(3);
  vector<double>    nA(3);
  double            nAlength = 0, nAn = 0, 
                    angle1 = 0, angle2 = 0, angle3 = 0;

  for(int i = 0; i < 3; i++)
    e0[i] = rawPointList_[rawTriangleList_[triangleId][1]][i]
      - rawPointList_[rawTriangleList_[triangleId][0]][i];

  for(int i = 0; i < 3; i++)
    e1[i] = (*p)[i] - rawPointList_[rawTriangleList_[triangleId][0]][i];

  // the point is the vertex 0
  if((!e1[0])&&(!e1[1])&&(!e1[2])) return true;
  
  // the point is the vertex 1
  if((e1[0] == e0[0])&&(e1[1] == e0[1])&&(e1[2] == e0[2])) return true;

  // cross product
  nAlength = 0;
  for(int i = 0; i < 3; i++){
    nA[i] = e0[(i+1)%3]*e1[(i+2)%3] - e0[(i+2)%3]*e1[(i+1)%3];
    nAlength += nA[i]*nA[i];
  }
  nAlength = sqrt(nAlength);

  // dot product
  angle1 = 0;
  for(int i = 0; i < 3; i++) angle1 += n[i]*nA[i];

  angle1 /= (nLength*nAlength);

  if(isnan(angle1)) return false;

  // now do exactly the same thing for the other vertices
  for(int i = 0; i < 3; i++)
    e0[i] = rawPointList_[rawTriangleList_[triangleId][2]][i]
      - rawPointList_[rawTriangleList_[triangleId][1]][i];

  for(int i = 0; i < 3; i++)
    e1[i] = (*p)[i] - rawPointList_[rawTriangleList_[triangleId][1]][i];

  // the point is the vertex 0
  if((!e1[0])&&(!e1[1])&&(!e1[2])) return true;
  
  // the point is the vertex 1
  if((e1[0] == e0[0])&&(e1[1] == e0[1])&&(e1[2] == e0[2])) return true;

  // cross product
  nAlength = 0;
  for(int i = 0; i < 3; i++){
    nA[i] = e0[(i+1)%3]*e1[(i+2)%3] - e0[(i+2)%3]*e1[(i+1)%3];
    nAlength += nA[i]*nA[i];
  }
  nAlength = sqrt(nAlength);

  // dot product
  angle2 = 0;
  for(int i = 0; i < 3; i++) angle2 += n[i]*nA[i];

  angle2 /= (nLength*nAlength);

  if(isnan(angle2)) return false;
 
  // last edge now

  for(int i = 0; i < 3; i++)
    e0[i] = rawPointList_[rawTriangleList_[triangleId][0]][i]
      - rawPointList_[rawTriangleList_[triangleId][2]][i];

  for(int i = 0; i < 3; i++)
    e1[i] = (*p)[i] - rawPointList_[rawTriangleList_[triangleId][2]][i];

  // the point is the vertex 0
  if((!e1[0])&&(!e1[1])&&(!e1[2])) return true;
  
  // the point is the vertex 1
  if((e1[0] == e0[0])&&(e1[1] == e0[1])&&(e1[2] == e0[2])) return true;

  // cross product
  nAlength = 0;
  for(int i = 0; i < 3; i++){
    nA[i] = e0[(i+1)%3]*e1[(i+2)%3] - e0[(i+2)%3]*e1[(i+1)%3];
    nAlength += nA[i]*nA[i];
  }
  nAlength = sqrt(nAlength);

  // dot product
  angle3 = 0;
  for(int i = 0; i < 3; i++) angle3 += n[i]*nA[i];

  angle3 /= (nLength*nAlength);

  if(isnan(angle3)) return false;



  if((angle1 > 0.66)&&(angle2 > 0.66)&&(angle3 > 0.66)){
    return true;
  }

  return false;
}

inline int TriangleMesh::projectPointOnTrianglePlane(const int &triangleId, 
  const vector<double> *point, vector<double> *projectedPoint) const{
      
  if((triangleId < 0)||(triangleId >= triangles_.size())) return -1;

  if((!point)||(point->size() != 3)) return -2;

  if(!projectedPoint) return -3;

  projectedPoint->resize(3);

  vector<double>  n(3);
  double          nLength = 0;
  getTriangleNormal(triangleId, &n);
  for(int i = 0; i < 3; i++) nLength += n[i]*n[i];
  nLength = sqrt(nLength);

  vector<double>  p0p(3);
  double          p0pLength = 0;

  for(int i = 0; i < 3; i++){
    p0p[i] = (*point)[i] - rawPointList_[rawTriangleList_[triangleId][0]][i];
    p0pLength += p0p[i]*p0p[i];
  }
  p0pLength = sqrt(p0pLength);

  double angle = 0;

  // dot product
  for(int i = 0; i < 3; i++) angle += p0p[i]*n[i];
  angle /= (p0pLength*nLength);

  if(isnan(angle)) return -4;

  double distance = angle*p0pLength;

  vector<double>  projection(3);
  for(int i = 0; i < 3; i++)
    projection[i] = -distance*n[i]/nLength;

  for(int i = 0; i < 3; i++)
    (*projectedPoint)[i] = (*point)[i] + projection[i];

  return 0;
}

inline int TriangleMesh::setVertexPoint(const int &vertexId, 
  const vector<double> &p){

  if((vertexId < 0)||(vertexId >= rawPointList_.size())) return -1;

  rawPointList_[vertexId] = p;

  return 0;
}

#endif
