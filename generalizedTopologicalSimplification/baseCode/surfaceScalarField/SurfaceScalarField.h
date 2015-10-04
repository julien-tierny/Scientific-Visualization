/*
 * file:                SurfaceScalarField.h
 * description:         Scalar fields on PL 2-manifold handling.
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

#ifndef                 _SCALAR_FIELD_SURFACE_H
#define                 _SCALAR_FIELD_SURFACE_H

#include                <Debug.h>
#include                <TriangleMesh.h>

#include                <algorithm>
#include                <set>
#include                <iomanip>

#include                <math.h>

using namespace std;

class SurfaceScalarField : public Debug{

  public:

    enum VertexType {Regular, Minimum, Saddle, Maximum};

    SurfaceScalarField();

    void operator<<(istream &f);

    void operator>>(ostream &f);

    inline  operator const vector<double>* () const{
      return (const vector<double> *) &vertexScalars_;};

    inline const vector<int>& getDebugExtraCriticalPoints() const{
      return debugExtraCriticalPoints_;};

    inline int getEdgeIsoContourTime(const int &triangleId, 
      const int &edgeTriangleId, const double &isoValue,
      double &contourTime) const;

    int getIsoContour(const int &seedVertexId, 
      vector<pair<pair<int, int>, double> > &contour) const;

    int getIsoContour(const int &seedTriangleId, const double &isoValue,
      vector<pair<pair<int, int>, double> > &contour) const;

    inline const vector<int>* getMaximumList(){
      if(needToUpdateCriticalPoints_)
        updateCriticalPoints();
      return (const vector<int> *) &maximumList_;
    }
    
    inline const double getMaximumValue() const { return maxValue_;};

    inline const vector<int>* getMinimumList(){
      if(needToUpdateCriticalPoints_)
        updateCriticalPoints();
      return (const vector<int> *) &minimumList_;
    }

    inline double getMinimumValue() const { return minValue_;};

    inline int getNumberOfScalars() const { return vertexScalars_.size();};

    int getPointScalar(const int &triangleId, 
      const vector<double> &baryCentrics, double &pointScalar) const{};

    inline const vector<int>* getSaddleList(){
      if(needToUpdateCriticalPoints_)
        updateCriticalPoints();
      return (const vector<int> *) &saddleList_;
    }

    inline const TriangleMesh* getTriangleMesh() const {return triangleMesh_;};

    inline int getVertexNormalizedScalar(const int &vertexId, 
      double &vertexNormalizedScalar) const;
    
    inline int getVertexScalar(const int &vertexId, double &vertexScalar) const;

    inline const vector<int>* getVertexTypes() {
      if(needToUpdateCriticalPoints_) updateCriticalPoints();
      return &vertexTypes_;};

    inline bool isEdgeOnIsoContour(const int &triangleId, 
      const int &triangleEdgeId, const double &isoValue) const;

    inline bool isSosLowerThan(const int &vertexId0, const int &vertexId1)const;

    inline bool isSosHigherThan(const int &vertexId0, 
      const int &vertexId1) const;

    inline bool isTriangleOnIsoContour(const int &triangleId, 
      const double &isoValue) const;

    int registerTriangleMesh(TriangleMesh *triangleMesh);

    inline int reset();

    int saveTopology(const string &path) const;

    inline int setVertexScalar(const int &vertexId, 
      const double &vertexScalar);

    int simplificationCheck(const vector<double> &originalField,
      const double &originalGlobalMin, const double &originalGlobalMax,
      const vector<int> &authorizedMinima, 
      const vector<int> &originalSaddleList,
      const vector<int> &authorizedMaxima);

    int simplify(const vector<int> &authorizedMinima, 
      const vector<int> &authorizedMaxima);

    int updateCriticalPoints();

  protected:

    bool                    needToUpdateCriticalPoints_;
    double                  minValue_, maxValue_;

    vector<int>             minimumList_,
                            saddleList_,
                            maximumList_;

    vector<int>             debugExtraCriticalPoints_;

    vector<int>             saddleMultiplicity_;

    vector<int>             vertexTypes_;

    TriangleMesh            *triangleMesh_;

    vector<double>          vertexScalars_;
    vector<int>             vertexSoSoffsets_;
};

inline int SurfaceScalarField::getEdgeIsoContourTime(const int &triangleId, 
  const int &edgeTriangleId, const double &isoValue, double &contourTime) const{

  if(!triangleMesh_) return -1;

  if(!vertexScalars_.size()) return -2;

  if((triangleId < 0)||(triangleId >= triangleMesh_->getNumberOfTriangles()))
    return -3;

  if((edgeTriangleId < 0)||(edgeTriangleId >= 3))
    return -4;

  const Triangle  *t = triangleMesh_->getTriangle(triangleId);
  const Edge      *e = t->getEdge(edgeTriangleId);
  int             v0, v1;

  e->getVertexId(0, v0);
  e->getVertexId(1, v1);

  double div = fabs(vertexScalars_[v1] - vertexScalars_[v0]);
  if(div < 0.000001) div = 0.000001;

  contourTime = fabs(isoValue - vertexScalars_[v0])/div;

  return 0;
}

inline int SurfaceScalarField::getVertexNormalizedScalar(const int &vertexId,
  double &vertexNormalizedScalar) const{
  
  if((vertexId < 0)||(vertexId >= vertexScalars_.size())){
    stringstream msg;
    msg << "[SurfaceScalarField] Vertex Id #" << vertexId
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return -1;
  }

  if((minValue_ == -1)&&(minValue_ == maxValue_))
    return -2;

  vertexNormalizedScalar = 
    (vertexScalars_[vertexId] - minValue_)/(maxValue_ - minValue_);
  
  return 0;  
}

inline int SurfaceScalarField::getVertexScalar(const int &vertexId,
  double &vertexScalar) const{
  
  if((vertexId < 0)||(vertexId >= vertexScalars_.size())){
    stringstream msg;
    msg << "[SurfaceScalarField] Vertex Id #" << vertexId
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return -1;
  }

  vertexScalar = vertexScalars_[vertexId];
  
  return 0;  
}

inline bool SurfaceScalarField::isEdgeOnIsoContour(const int &triangleId, 
  const int &triangleEdgeId, const double &isoValue) const{

  if(!triangleMesh_) return false;

  if(!vertexScalars_.size()) return false;

  if((triangleId < 0)||(triangleId >= triangleMesh_->getNumberOfTriangles()))
    return false;

  if((triangleEdgeId < 0)||(triangleEdgeId >= 3)) return false;

  const Triangle  *t = triangleMesh_->getTriangle(triangleId);
  const Edge      *e = t->getEdge(triangleEdgeId);
  int             v0, v1;

  e->getVertexId(0, v0);
  e->getVertexId(1, v1);

  return (((vertexScalars_[v0] < isoValue)&&(vertexScalars_[v1] >= isoValue))
    ||((vertexScalars_[v1] < isoValue)&&(vertexScalars_[v0] >= isoValue)));

  return 0;  
}

inline bool SurfaceScalarField::isSosHigherThan(const int &vertexId0, 
  const int &vertexId1)const{
 
  return ((vertexScalars_[vertexId0] > vertexScalars_[vertexId1])
    ||((vertexScalars_[vertexId0] == vertexScalars_[vertexId1])
      &&(vertexSoSoffsets_[vertexId0] > vertexSoSoffsets_[vertexId1])));
}

inline bool SurfaceScalarField::isSosLowerThan(const int &vertexId0, 
  const int &vertexId1)const{
  
  return ((vertexScalars_[vertexId0] < vertexScalars_[vertexId1])
    ||((vertexScalars_[vertexId0] == vertexScalars_[vertexId1])
      &&(vertexSoSoffsets_[vertexId0] < vertexSoSoffsets_[vertexId1])));
}

inline bool SurfaceScalarField::isTriangleOnIsoContour(const int &triangleId, 
  const double &isoValue) const{
      
  if(!triangleMesh_) return false;
  if((triangleId < 0)
    ||(triangleId >= triangleMesh_->getNumberOfTriangles())) return false;
  
  return ((isEdgeOnIsoContour(triangleId, 0, isoValue))
    ||(isEdgeOnIsoContour(triangleId, 1, isoValue))
    ||(isEdgeOnIsoContour(triangleId, 2, isoValue)));
}

inline int SurfaceScalarField::reset(){

  needToUpdateCriticalPoints_ = true;
  minValue_ = maxValue_ = -1;

  for(int i = 0; i < vertexSoSoffsets_.size(); i++)
    vertexSoSoffsets_[i] = i;

  return 0;
}

inline int SurfaceScalarField::setVertexScalar(const int &vertexId, 
  const double &vertexScalar){
  
  if((vertexId < 0)||(vertexId >= vertexScalars_.size())){
    stringstream msg;
    msg << "[SurfaceScalarField] Vertex Id #" << vertexId
      << " out of range." << endl;
    dMsg(cerr, msg.str(), 10);
    return -1;
  }

  vertexScalars_[vertexId] = vertexScalar;

  if((minValue_ == -1)&&(minValue_ == maxValue_)){
    minValue_ = vertexScalar;
    maxValue_ = vertexScalar;
  }
  else{
    if(vertexScalar < minValue_) minValue_ = vertexScalar;
    if(vertexScalar > maxValue_) maxValue_ = vertexScalar;
  }

  needToUpdateCriticalPoints_ = true;

  return 0;  
}

#endif
