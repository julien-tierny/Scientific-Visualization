/*
 * file:                SurfaceScalarField.cpp
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

#include                "SurfaceScalarField.h"

bool                    increasingOrder_ = true;

struct sweepCmp{
  bool operator() (const pair<double, pair<int, int> > &v0,
    const pair<double, pair<int, int> > &v1){

    if(increasingOrder_){
      return ((v0.first < v1.first)
        ||((v0.first == v1.first)
          &&(v0.second.first < v1.second.first)));
    }
    else{
      return ((v0.first > v1.first)
        ||((v0.first == v1.first)
          &&(v0.second.first > v1.second.first)));
    }
  };
} sortCmp;

SurfaceScalarField::SurfaceScalarField(){

  triangleMesh_ = NULL;

  minValue_ = -1; 
  maxValue_ = -1;

  needToUpdateCriticalPoints_ = false;
}

void SurfaceScalarField::operator<<(istream &f){

  int     vertexNumber = 0;
  double  vertexScalar = 0;
  string keyword;

  if(!triangleMesh_)
    dMsg(cout, 
  "[SurfaceScalarField] Warning, reading scalar data but no mesh registered!\n",
    1);

  f >> keyword;
  if(keyword != "SSF")
    dMsg(cerr, 
    "[SurfaceScalarField] Input file does not seem to be a valid *ssf file.\n",
    0);

  f >> vertexNumber;

  if((triangleMesh_)&&(vertexNumber != triangleMesh_->getNumberOfVertices())){
    stringstream msg;
    msg << "[SurfaceScalarField] Warning! Reading " << vertexNumber 
      << " scalars but the mesh has " << triangleMesh_->getNumberOfVertices() 
      << " vertices!" << endl;
    dMsg(cerr, msg.str(), 0);
  }

  vertexScalars_.resize(vertexNumber);
  vertexTypes_.resize(vertexNumber);

  for(int i = 0; i < vertexNumber; i++){
    // there's no hope to get a decent enough precision with doubles if we're 
    // asking for a 15 digit precision. for safety, let's ask 10.
    f >> setprecision(10) >> vertexScalar;
    setVertexScalar(i, vertexScalar);
  }
}

void SurfaceScalarField::operator>>(ostream &f){

  f << "SSF" << endl;
  f << vertexScalars_.size() << endl;

  // transform the symbolic perturbation into a numerical perturbation
  vector<pair<double, pair<int, int> > > perturbation(vertexScalars_.size());
  for(int i = 0; i < vertexScalars_.size(); i++){
    perturbation[i].first = vertexScalars_[i];
    perturbation[i].second.first = vertexSoSoffsets_[i];
    perturbation[i].second.second = i;
  }

  // sort the thing
  increasingOrder_ = true;  
  sort(perturbation.begin(), perturbation.end(), sortCmp);

  // transform the perturbation 
  for(int i = 0; i < perturbation.size(); i++){
    if(i){
      if(perturbation[i].first <= perturbation[i - 1].first){
        perturbation[i].first = perturbation[i - 1].first
          // 10 digit precision.
          //
          // there's no hope to get a decent enough precision with doubles 
          // if we're asking for a 15 digit precision (theoretical bound). 
          // for safety, let's ask 10.
          //
          // with values spanning from 0 to 1, it'll work until we have more
          // than 10 power of 10 vertices in there (several billions, I don't
          // think so :) ).
          + 0.000000001;
      }
    }
    vertexScalars_[perturbation[i].second.second] = perturbation[i].first;

    // update for normalization purpose
    if((!i)||(vertexScalars_[perturbation[i].second.second] > maxValue_))
      maxValue_ = vertexScalars_[perturbation[i].second.second];
    
    if((!i)||(vertexScalars_[perturbation[i].second.second] < minValue_))
      minValue_ = vertexScalars_[perturbation[i].second.second];
  }

  for(int i = 0; i < vertexScalars_.size(); i++)
    f << setprecision(10) << vertexScalars_[i] << endl;
}

int SurfaceScalarField::getIsoContour(const int &seedVertexId, 
  vector<pair<pair<int, int>, double> > &contour) const{

  if(!triangleMesh_) return -1;

  if(!vertexScalars_.size()) return -2;

  if((seedVertexId < 0)||(seedVertexId >= vertexScalars_.size())) return -3;

  int             tId;
  const Vertex    *v = triangleMesh_->getVertex(seedVertexId);
  const Triangle  *t;

  for(int i = 0; i < v->getNumberOfTriangles(); i++){
    // find a triangle on the contour
    v->getTriangleId(i, tId);
    if(isTriangleOnIsoContour(tId, vertexScalars_[seedVertexId]))
      return getIsoContour(tId, vertexScalars_[seedVertexId], contour);
  }

  return 0;
}

int SurfaceScalarField::getIsoContour(const int &triangleId, 
  const double &isoValue, vector<pair<pair<int, int>, double> > &contour) const{

  DebugTimer timer;

  if(!triangleMesh_) return -1;

  if(!vertexScalars_.size()) return -2;

  if((triangleId < 0)||
    (triangleId >= triangleMesh_->getNumberOfTriangles())) return -3;

  vector<bool> visitedTriangles(triangleMesh_->getNumberOfTriangles());
  for(int i = 0; i < visitedTriangles.size(); i++) visitedTriangles[i] = false;

  const Edge                    *e = NULL;
  const Triangle                *t = NULL;
  pair<pair<int, int>, double>  contourPoint;
  int                           nextTriangleId, edgeId, pivotEdgeId = -1;
  
  contour.clear();

  nextTriangleId = triangleId;

  do{
    for(int i = 0; i < 3; i++){
      t = triangleMesh_->getTriangle(nextTriangleId);
      t->getEdgeId(i, edgeId);
      if((isEdgeOnIsoContour(nextTriangleId, i, isoValue))
        &&(edgeId != pivotEdgeId)){
        contourPoint.first.first = nextTriangleId;
        contourPoint.first.second = i;
        getEdgeIsoContourTime(nextTriangleId, i, isoValue, contourPoint.second);
        
        contour.push_back(contourPoint);
        visitedTriangles[nextTriangleId] = true;
      
        // put also the other edge of the first triangle
        if(nextTriangleId == triangleId){
          for(int j = 0; j < 3; j++){
            if((j != i)&&(isEdgeOnIsoContour(triangleId, j, isoValue))){
              contourPoint.first.first = triangleId;
              contourPoint.first.second = j;
              getEdgeIsoContourTime(triangleId, j, isoValue, 
                contourPoint.second);
              contour.push_back(contourPoint);
              break;
            }
          }
        }

        // get the next triangle to visit
        t = triangleMesh_->getTriangle(contourPoint.first.first);
        t->getEdgeId(contourPoint.first.second, pivotEdgeId);
        e = t->getEdge(contourPoint.first.second);
        e->getTriangleId(0, nextTriangleId);
        if(nextTriangleId == contourPoint.first.first)
          e->getTriangleId(1, nextTriangleId);
        break;
      }
    }
  }while((nextTriangleId != -1)&&(!visitedTriangles[nextTriangleId]));
 
  if(nextTriangleId == -1){


    vector<pair<pair<int, int>, double> > otherSegment;
    // we hit a boundary, we need to restart from the startTriangle...
    // then revert other segment and concat the two...
  
    nextTriangleId = triangleId;

    // retrieve the starting edge 
    t = triangleMesh_->getTriangle(nextTriangleId);
    t->getEdgeId(contour[0].first.second, pivotEdgeId);
    e = t->getEdge(contour[0].first.second);
   
    // find the next triangle in the other direction
    e->getTriangleId(0, nextTriangleId);
    if(nextTriangleId == contour[0].first.first)
      e->getTriangleId(1, nextTriangleId);


    if(nextTriangleId != -1){
      // and restart the process
      do{
        for(int i = 0; i < 3; i++){
          t = triangleMesh_->getTriangle(nextTriangleId);
          t->getEdgeId(i, edgeId);
          if((isEdgeOnIsoContour(nextTriangleId, i, isoValue))
            &&(edgeId != pivotEdgeId)){
              contourPoint.first.first = nextTriangleId;
            contourPoint.first.second = i;
            getEdgeIsoContourTime(nextTriangleId, i, isoValue, 
              contourPoint.second);
            otherSegment.push_back(contourPoint);
            visitedTriangles[nextTriangleId] = true;
      
            // get the next triangle to visit
            t = triangleMesh_->getTriangle(contourPoint.first.first);
            t->getEdgeId(contourPoint.first.second, pivotEdgeId);
            e = t->getEdge(contourPoint.first.second);
            e->getTriangleId(0, nextTriangleId);
            if(nextTriangleId == contourPoint.first.first)
              e->getTriangleId(1, nextTriangleId);
            break;
          }
        }
      }while((nextTriangleId != -1)&&(!visitedTriangles[nextTriangleId]));

      // now let's fill the final structure
      vector<pair<pair<int, int>, double> > tmpContour = contour;

      contour.resize(otherSegment.size() + tmpContour.size());
    
      int contourIndex = 0;
      for(int i = otherSegment.size() - 1; i >= 0; i--){
        contour[contourIndex] = otherSegment[i];
        contourIndex++;
      }
      for(int i = 0; i < tmpContour.size(); i++){
        contour[contourIndex] = tmpContour[i];
        contourIndex++;
      }  
    }
  }

  stringstream msg;
  msg << "[SurfaceScalarField] Contour extracted (" << contour.size()
    << " triangles) in " 
    << timer.getElapsedTime()
    << " s." << endl;

  dMsg(cout, msg.str(), 3);

  return 0;
}

int SurfaceScalarField::registerTriangleMesh(TriangleMesh *triangleMesh){

  if(!triangleMesh) return -1;

  triangleMesh_ = triangleMesh;

  vertexScalars_.resize(triangleMesh_->getNumberOfVertices());
  vertexSoSoffsets_.resize(triangleMesh_->getNumberOfVertices());
  // init offset
  for(int i = 0; i < vertexSoSoffsets_.size(); i++)
    vertexSoSoffsets_[i] = i;
  vertexTypes_.resize(triangleMesh_->getNumberOfVertices());

  return 0;
}

int SurfaceScalarField::saveTopology(const string &path) const{

  ofstream f(path.data(), ios::out);

  if(!f){
    stringstream msg;
    msg << "[SurfaceScalarField] Cannot write in file '" 
      << path << "'." << endl;
    dMsg(cerr, msg.str(), 1);
    return -1;
  }

  f << "STF" << endl;

  f << minimumList_.size() << " " 
    << maximumList_.size() << endl;

  for(int i = 0; i < minimumList_.size(); i++){
    f << "m " << minimumList_[i] << endl;
  }

  for(int i = 0; i < maximumList_.size(); i++){
    f << "M " << maximumList_[i] << endl;
  }

  f.close();

  return 0;
}

int SurfaceScalarField::simplificationCheck(const vector<double> &originalField,
  const double &originalGlobalMin, const double &originalGlobalMax,
  const vector<int> &authorizedMinima, const vector<int> &originalSaddleList,
  const vector<int> &authorizedMaxima){

  double maxDistance = 0, l1Norm = 0, l2Norm = 0;
  double normalizedValue = 0, originalNormalizedValue = 0;

    
  // get the actual delta of the simplification
  stringstream deltaMsg;

  for(int i = 0; i < vertexScalars_.size(); i++){
    if((!i)||(fabs(vertexScalars_[i] - originalField[i]) > maxDistance))
      maxDistance = fabs(vertexScalars_[i] - originalField[i]);
    l1Norm += fabs(vertexScalars_[i] - originalField[i]);
    l2Norm += (vertexScalars_[i] - originalField[i])
      *(vertexScalars_[i] - originalField[i]);
  }
  l2Norm = sqrt(l2Norm);

  deltaMsg
    << "[SurfaceScalarField] Simplification L_1 norm: "
    << l1Norm
    << " ("
    << (l1Norm/vertexScalars_.size())*100
      /(originalGlobalMax - originalGlobalMin)
    << " %)"
    << endl;

  deltaMsg
    << "[SurfaceScalarField] Simplification L_2 norm: "
    << l2Norm
    << endl;

  deltaMsg 
    << "[SurfaceScalarField] Simplification L_infinity norm: " 
    << maxDistance 
    << " ("
    << maxDistance*100/(originalGlobalMax - originalGlobalMin)
    << " %)"
    << endl;
  dMsg(cout, deltaMsg.str(), 2);

  // make sure all of the new critical points are authorized.
  stringstream criticalMsg;

  debugExtraCriticalPoints_.clear();
  updateCriticalPoints();
    
  bool found = false, criticalOk = true;
  for(int i = 0; i < minimumList_.size(); i++){
    found = false;
    for(int j = 0; j < authorizedMinima.size(); j++){
      if(minimumList_[i] == authorizedMinima[j]){
        found = true;
        if(originalField[minimumList_[i]] != vertexScalars_[minimumList_[i]]){
          stringstream msg;
          msg << "[SurfaceScalarField] ";
          if(triangleMesh_->getVertex(minimumList_[i])->isOnBoundary())
            msg << "Boundary min";
          else 
            msg << "Min";
          msg << " " << minimumList_[i] 
            << " shifted from " << originalField[minimumList_[i]] 
            << " to " << vertexScalars_[minimumList_[i]] << endl;
          dMsg(cout, msg.str(), 3);
        }
        break;
      }
    }
    if(!found){
      criticalOk = false;
      if(debugLevel_ >= 2){
        criticalMsg << "[SurfaceScalarField] Min vertex " << 
          minimumList_[i] << " not authorized!" << endl;
      }
      debugExtraCriticalPoints_.push_back(minimumList_[i]);
    }
  }
  for(int i = 0; i < saddleList_.size(); i++){
    found = false;
    for(int j = 0; j < originalSaddleList.size(); j++){
      if(saddleList_[i] == originalSaddleList[j]){
        found = true;
        if(originalField[saddleList_[i]] != vertexScalars_[saddleList_[i]]){
          stringstream msg;
          msg << "[SurfaceScalarField] ";
          if(triangleMesh_->getVertex(saddleList_[i])->isOnBoundary())
            msg << "Boundary saddle";
          else
            msg << "Saddle";
          msg << " " << saddleList_[i] 
            << " shifted from " << originalField[saddleList_[i]] 
            << " to " << vertexScalars_[saddleList_[i]] << endl;
          dMsg(cout, msg.str(), 3);
        }
        break;
      }
    }
    if(!found){
      if(debugLevel_ >= 3){
        criticalMsg << "[SurfaceScalarField] ";
        if(triangleMesh_->getVertex(saddleList_[i])->isOnBoundary()){
          criticalMsg << "Boundary saddle";
        }
        else criticalMsg << "Saddle";
        criticalMsg << " emerged in vertex " << 
          saddleList_[i] << "." << endl;
      }
    }
  }
  for(int i = 0; i < maximumList_.size(); i++){
    found = false;
    for(int j = 0; j < authorizedMaxima.size(); j++){
      if(maximumList_[i] == authorizedMaxima[j]){
        found = true;
        if(originalField[maximumList_[i]] != vertexScalars_[maximumList_[i]]){
          stringstream msg;
          msg << "[SurfaceScalarField] ";
          if(triangleMesh_->getVertex(maximumList_[i])->isOnBoundary())
            msg << "Boundary max";
          else
            msg << "Max";
          msg << " " << maximumList_[i] 
            << " shifted from " << originalField[maximumList_[i]] 
            << " to " << vertexScalars_[maximumList_[i]] << endl;
          dMsg(cout, msg.str(), 3);
        }
        break;
      }
    }
    if(!found){
      criticalOk = false;
      if(debugLevel_ >= 2){
        criticalMsg << "[SurfaceScalarField] Max vertex " << 
          maximumList_[i] << " not authorized!" << endl;
      }
      debugExtraCriticalPoints_.push_back(maximumList_[i]);
    }
  }

  if((minimumList_.size() > authorizedMinima.size())
    ||(maximumList_.size() > authorizedMaxima.size()))
    criticalOk = false;
  else if((minimumList_.size() < authorizedMinima.size())
    ||(maximumList_.size() < authorizedMaxima.size())){
    stringstream msg;
    msg << "[SurfaceScalarField] Constraints deleted. m: " 
      << minimumList_.size() << "/" << authorizedMinima.size() << ", M: "
      << maximumList_.size() << "/" << authorizedMaxima.size()
      << endl;
    dMsg(cout, msg.str(), 2);
  }

  if(!criticalOk){
    criticalMsg 
<< "[SurfaceScalarField] Field simplification inserted new critical points..." 
      << endl;
  }
  else{
    criticalMsg <<
      "[SurfaceScalarField] Constrained topology correctly enforced." << endl;
  }

  dMsg(cout, criticalMsg.str(), 2);
}

int SurfaceScalarField::simplify(const vector<int> &authorizedMinima,
  const vector<int> &authorizedMaxima){

  DebugTimer timer;

  if((!authorizedMinima.size())||(!authorizedMaxima.size()))
    return -1;

  int             iterationNb = 0;
  double          originalGlobalMin, originalGlobalMax;
  vector<double>  originalField;
  vector<int>     originalMinimumList, originalSaddleList, originalMaximumList;
  vector<bool>    authorizedExtrema(
                    triangleMesh_->getNumberOfVertices(), false);
 
  // stats purpose
  if(needToUpdateCriticalPoints_) updateCriticalPoints();

  originalMinimumList = minimumList_;
  originalSaddleList = saddleList_;
  originalMaximumList = maximumList_;
  //

  if(debugLevel_ >= 2){
    originalField = vertexScalars_;
    originalGlobalMin = minValue_;
    originalGlobalMax = maxValue_;
  }

  bool needForMoreIterations = true;

  do{

    // one pass for the sub-level sets, one for the sur-level sets
    for(int j = 0; j < 2; j++){

      if(!j) increasingOrder_ = true;
      else increasingOrder_ = false;

      // should be allocated to false
      vector<bool>            visitedVertices(
                                triangleMesh_->getNumberOfVertices(), false);

      // first:  scalarValue
      // second.first: SoS offset
      // second.second: vertex index
      set<pair<double, pair<int, int> >, sweepCmp > sweepFront;

      int                     adjustmentPos = 0;
      vector<int>             adjustmentSequence(
                                triangleMesh_->getNumberOfVertices());

      int     vertexId = -1, nId = -1, offset = 0;
      const   Vertex *v = NULL;


      // 1) re-construct the topology of the sub level sets
      // add the authorized minima to the front
      if(increasingOrder_){
        for(int i = 0; i < authorizedMinima.size(); i++){
          authorizedExtrema[authorizedMinima[i]] = true;
          sweepFront.insert(
            pair<double, pair<int, int> >(
              vertexScalars_[authorizedMinima[i]],
              pair<int, int>(vertexSoSoffsets_[authorizedMinima[i]],
                authorizedMinima[i])));
          visitedVertices[authorizedMinima[i]] = true;
        }
      }
      else{
        for(int i = 0; i < authorizedMaxima.size(); i++){
          authorizedExtrema[authorizedMaxima[i]] = true;
          sweepFront.insert(
            pair<double, pair<int, int> >(
              vertexScalars_[authorizedMaxima[i]],
              pair<int, int>(vertexSoSoffsets_[authorizedMaxima[i]],
                authorizedMaxima[i])));
          visitedVertices[authorizedMaxima[i]] = true;
        }
      }

      // constrained sweep loop
      // ---------------------------
      //  
      // sweepFront.begin() will always point to the lowest vertex in the
      // front.

      do{
        // N steps (N: vertices)

        vertexId = sweepFront.begin()->second.second;
        sweepFront.erase(sweepFront.begin());

        v = triangleMesh_->getVertex(vertexId);
        if(!v) return -2;

        // this loop can be considered to take a constant number of steps
        // (the number of vertices in the star is always many orders of 
        // magnitude lower than the total number of vertices)
        for(int i = 0; i < v->getNumberOfStarVertices(); i++){

          v->getStarVertexId(i, nId);
          if(nId == -1) return -3;

          if(!visitedVertices[nId]){
            // log(n) (n: vertices in the front, n << N)
            sweepFront.insert(
              pair<double, pair<int, int> >(
                vertexScalars_[nId],
                  pair<int, int>(vertexSoSoffsets_[nId], nId)));

            visitedVertices[nId] = true;
          }
        }

        adjustmentSequence[adjustmentPos] = vertexId;
        adjustmentPos++;

      }while(!sweepFront.empty());

      // now enforce the monocity of the field with regard to the constrained
      // sweep loop
      
      if(increasingOrder_)
        offset = 0;
      else
        offset = adjustmentSequence.size() + 1; 

      for(int i = 0; i < adjustmentSequence.size(); i++){

        if(increasingOrder_){

          if((i)
            &&(vertexScalars_[adjustmentSequence[i]] 
              <= vertexScalars_[adjustmentSequence[i - 1]])){
 
            vertexScalars_[adjustmentSequence[i]] = 
              vertexScalars_[adjustmentSequence[i - 1]];
          }
          offset++;
        }
        else{
          if((i)
            &&(vertexScalars_[adjustmentSequence[i]] 
              >= vertexScalars_[adjustmentSequence[i - 1]])){
    
            vertexScalars_[adjustmentSequence[i]]
              = vertexScalars_[adjustmentSequence[i - 1]];
          }
          offset--;
        }
 
        vertexSoSoffsets_[adjustmentSequence[i]] = offset;

        // update global extrema for normalization purpose
        if((!i)||(vertexScalars_[adjustmentSequence[i]] < minValue_)){
          minValue_ = vertexScalars_[adjustmentSequence[i]];
        }
        if((!i)||(vertexScalars_[adjustmentSequence[i]] > maxValue_)){
          maxValue_ = vertexScalars_[adjustmentSequence[i]];
        }

      }
    }

    iterationNb++;

    // now check if we need an extra pass
    updateCriticalPoints();

    needForMoreIterations = false;

    if(maximumList_.size() > authorizedMaxima.size()) 
      needForMoreIterations = true;
    if(minimumList_.size() > authorizedMinima.size()) 
      needForMoreIterations = true;

    if(!needForMoreIterations){

      for(int i = 0; i < minimumList_.size(); i++){
        if(!authorizedExtrema[minimumList_[i]]){
          needForMoreIterations = true;
          break;
        }
      }
      if(!needForMoreIterations){
        for(int i = 0; i < maximumList_.size(); i++){
          if(!authorizedExtrema[maximumList_[i]]){
            needForMoreIterations = true;
            break;
          }
        }
      }
    }

  }while(needForMoreIterations);

  stringstream msg;

  msg << "[SurfaceScalarField] Field simplified (ext.: "
    << originalMinimumList.size() + originalMaximumList.size()
    << " -> "
    << authorizedMinima.size() + authorizedMaxima.size()
    << ") in _ "
    << timer.getElapsedTime()
    << " _ s. (" << iterationNb << " it)"<< endl;

  dMsg(cout, msg.str(), 1);

  if(debugLevel_ >= 2)
    simplificationCheck(originalField, originalGlobalMin, originalGlobalMax,
      authorizedMinima, originalSaddleList, authorizedMaxima);

  return 0;
}

int SurfaceScalarField::updateCriticalPoints(){

  if(!triangleMesh_) return -1;

  if(!vertexScalars_.size()) return -2;

  DebugTimer timer;

  int           nId, signChange, extraMultiplicity = 0, multipleSaddles = 0;
  bool          isMinimum, isMaximum, prevIsLower, firstIsLower;
  const Vertex  *v = NULL;

  minimumList_.clear();
  maximumList_.clear();
  saddleList_.clear();
  saddleMultiplicity_.clear();

  for(int i = 0; i < triangleMesh_->getNumberOfVertices(); i++){

    isMinimum = true;
    isMaximum = true;
    signChange = 0;

    v = triangleMesh_->getVertex(i);

    for(int j = 0; j < v->getNumberOfStarVertices(); j++){
      v->getStarVertexId(j, nId);

      if(isSosLowerThan(i, nId)){
        isMaximum = false;

        if((j)&&(prevIsLower)) signChange++;

        if(!j) firstIsLower = false;

        if((j == v->getNumberOfStarVertices() - 1)&&(firstIsLower)
          &&(!v->isOnBoundary())){
          // is there a sign change at the end of the loop.
          signChange++;
        }
          
        prevIsLower = false;
      }
      if(isSosHigherThan(i, nId)){
        isMinimum = false;

        if((j)&&(!prevIsLower)) signChange++;

        if(!j) firstIsLower = true;
        
        if((j == v->getNumberOfStarVertices() - 1)&&(!firstIsLower)
          &&(!v->isOnBoundary())){
          // is there a sign change at the end of the loop.
          signChange++;
        }

        prevIsLower = true;
      }
    }

    if(isMinimum){
      vertexTypes_[i] = Minimum;
      minimumList_.push_back(i);
      if(debugLevel_ >= 2){
        stringstream msg;
        msg << "[SurfaceScalarField] Vertex " << i << " is a min ("
          << vertexScalars_[i] << ", " << vertexSoSoffsets_[i] 
          << ")." << endl;
        dMsg(cout, msg.str(), 10);
      }
    }
    else if(isMaximum){
      vertexTypes_[i] = Maximum;
      maximumList_.push_back(i);
      if(debugLevel_ >= 2){
        stringstream msg;
        msg << "[SurfaceScalarField] Vertex " << i << " is a max ("
          << vertexScalars_[i] << ", " << vertexSoSoffsets_[i] 
          << ")." << endl;
        dMsg(cout, msg.str(), 10);
      }

    }
    else{
      if((v->isOnBoundary())&&(signChange == 1))
        vertexTypes_[i] = Regular;
      else if((!v->isOnBoundary())&&(signChange == 2))
        vertexTypes_[i] = Regular;
      else if(!v->isOnBoundary()){
        saddleList_.push_back(i);
        vertexTypes_[i] = Saddle;
        saddleMultiplicity_.push_back(signChange/2 - 1);
        extraMultiplicity += signChange/2 - 2;
        if(signChange/2 - 2) multipleSaddles++;

        if(debugLevel_ >= 2){
          stringstream msg;
          msg << "[SurfaceScalarField] Vertex " << i << " is a saddle ("
            << vertexScalars_[i] << ", " << vertexSoSoffsets_[i] 
            << ")." << endl;
          dMsg(cout, msg.str(), 10);
        }

      }
      else if(v->isOnBoundary()){
        saddleList_.push_back(i);
        vertexTypes_[i] = Saddle;
        saddleMultiplicity_.push_back(signChange - 1);
        extraMultiplicity += signChange - 2;
        if(signChange - 2) multipleSaddles++;

        if(debugLevel_ >= 2){
          stringstream msg;
          msg << "[SurfaceScalarField] Vertex " << i << " is a saddle ("
            << vertexScalars_[i] << ", " << vertexSoSoffsets_[i] 
            << ")." << endl;
          dMsg(cout, msg.str(), 10);
        }
      }
    }
  }

  stringstream msg;
  msg << "[SurfaceScalarField] Critical point update: "
    << timer.getElapsedTime()
    << " s. (" << minimumList_.size() << " m, "
    << saddleList_.size() << " s, "
    << maximumList_.size() << " M)" << endl;

  if((extraMultiplicity)&&(debugLevel_ >= 2)){
    msg << "[SurfaceScalarField] " << multipleSaddles
      << " multiple saddle(s) (total extra multiplicity " << extraMultiplicity
      << ")." << endl;
  }

  dMsg(cout, msg.str(), 1);

  // sanity check
  int chi = triangleMesh_->getNumberOfVertices() 
    - triangleMesh_->getNumberOfEdges()
    + triangleMesh_->getNumberOfTriangles();

  if(minimumList_.size() - (saddleList_.size() + extraMultiplicity) 
    + maximumList_.size() != chi){

    stringstream msg;
    msg << 
      "[SurfaceScalarField]" << 
      " Error, the Morse-Euler relation does not hold! (chi=" 
      << chi << ")." << endl;
    dMsg(cerr, msg.str(), 5);
  }

  needToUpdateCriticalPoints_ = false;

  return 0;
}
