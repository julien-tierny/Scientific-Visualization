/*
 * file:                Editor.cpp
 * description:         Class to handle data-structures and operations.
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

#include                "Editor.h"

int Editor::exportToVRML(const string &path){

  stringstream startMsg;
  startMsg << "[Editor] Saving scene in '" << path << "'..." << endl;
  dMsg(cout, startMsg.str(), 1);

  ofstream f(path.data(), ios::out);

  if(!f){
    stringstream msg;
    msg << "[Editor] Could not write file '" << path << "'." << endl;
    dMsg(cerr, msg.str(), 1);
    return -1;
  }

  f << "#VRML V2.0 utf8" << endl;
  f << "# File generated with `generalSimplification'." << endl;
  f << "# Copyright (C) Julien Tierny <tierny@telecom-paristech.fr> March 2012."
    << endl;

  f << "Group{" << endl;
  f << "\tchildren[" << endl;

  // feed the texture mapping with the scalar values
  vector<pair<double, double> > uvList(scalarField_.getNumberOfScalars());
  for(int i = 0; i < scalarField_.getNumberOfScalars(); i++){
    double vertexScalar = 0;
    scalarField_.getVertexScalar(i, vertexScalar);
    uvList[i] = pair<double, double>(
      (vertexScalar - initialMinValue_)/(initialMaxValue_ - initialMinValue_),
      0);
  }

  triangleMesh_.exportToVRML(f, uvList, "textures/texture.jpg");
  
  const vector<int> *minimumList = scalarField_.getMinimumList();
  for(int i = 0; i < minimumList->size(); i++){

    f << "# Minimum " << i << endl;
    f << "Transform{" << endl;
    f << "\ttranslation ";
    
    const vector<double> *p = 
      triangleMesh_.getVertex((*minimumList)[i])->getPoint();
    for(int j = 0; j < 3; j++) f << (*p)[j] << " ";
    f << endl;

    f << "\tchildren[" << endl;
    f << "\t\tShape{" << endl;
    f << "\t\t\tappearance Appearance {" << endl;
    f << "\t\t\t\tmaterial Material {" << endl;
    f << "\t\t\t\t\tdiffuseColor 0.254902 0.364706 0.690196" << endl;
    f << "\t\t\t\t\tspecularColor 0.8 0.8 0.8" << endl;
    f << "\t\t\t\t\tshininess 10" << endl;
    f << "\t\t\t\t}" << endl;
    f << "\t\t\t}" << endl;

    f << "\t\t\tgeometry Sphere { radius 0.05 }" << endl;

    f << "\t\t}" << endl;
    f << "\t]" << endl;
    f << "}" << endl;
  }

  const vector<int> *saddleList = scalarField_.getSaddleList();
  for(int i = 0; i < saddleList->size(); i++){
    
    f << "# Saddle " << i << endl;
    f << "Transform{" << endl;
    f << "\ttranslation ";
    
    const vector<double> *p = 
      triangleMesh_.getVertex((*saddleList)[i])->getPoint();
    for(int j = 0; j < 3; j++) f << (*p)[j] << " ";
    f << endl;

    f << "\tchildren[" << endl;
    f << "\t\tShape{" << endl;
    f << "\t\t\tappearance Appearance {" << endl;
    f << "\t\t\t\tmaterial Material {" << endl;
    f << "\t\t\t\t\tdiffuseColor 0.678431 0.701961 0.580392" << endl;
    f << "\t\t\t\t\tspecularColor 0.8 0.8 0.8" << endl;
    f << "\t\t\t\t\tshininess 10" << endl;
    f << "\t\t\t\t}" << endl;
    f << "\t\t\t}" << endl;

    f << "\t\t\tgeometry Sphere { radius 0.05 }" << endl;

    f << "\t\t}" << endl;
    f << "\t]" << endl;
    f << "}" << endl;
  }

  const vector<int> *maximumList = scalarField_.getMaximumList();
  for(int i = 0; i < maximumList->size(); i++){

    f << "# Maximum " << i << endl;
    f << "Transform{" << endl;
    f << "\ttranslation ";
    
    const vector<double> *p = 
      triangleMesh_.getVertex((*maximumList)[i])->getPoint();
    for(int j = 0; j < 3; j++) f << (*p)[j] << " ";
    f << endl;

    f << "\tchildren[" << endl;
    f << "\t\tShape{" << endl;
    f << "\t\t\tappearance Appearance {" << endl;
    f << "\t\t\t\tmaterial Material {" << endl;
    f << "\t\t\t\t\tdiffuseColor 0 0.321569 0" << endl;
    f << "\t\t\t\t\tspecularColor 0.8 0.8 0.8" << endl;
    f << "\t\t\t\t\tshininess 10" << endl;
    f << "\t\t\t\t}" << endl;
    f << "\t\t\t}" << endl;

    f << "\t\t\tgeometry Sphere { radius 0.05 }" << endl;

    f << "\t\t}" << endl;
    f << "\t]" << endl;
    f << "}" << endl;
  }

  f << "\t]" << endl;
  f << "}" << endl;  

  f.close();

  stringstream msg;
  msg << "[Editor] Scene saved in '" << path << "'." << endl;
  dMsg(cout, msg.str(), 1);

  return 0;
}

int Editor::loadScalarField(const string &scalarPath, 
  const string &topologyPath){

  if(!scalarPath.length()){
    stringstream msg;
    msg 
      << "[Editor] No user scalar field specified!" << endl;
    dMsg(cerr, msg.str(), 1);
    exit(-1);
  }
  else{

    ifstream f(scalarPath.data(), ios::in);

    if(!f){
      stringstream msg;
      msg << "[Editor] Cannot read file '" << scalarPath << "'!" << endl;
      dMsg(cerr, msg.str(), 1);
      exit(-1);
    }
    else{

      scalarField_ << f;

      stringstream msg;
      msg << "[Editor] Scalar field '" << scalarPath << "' read." << endl;
      dMsg(cout, msg.str(), 1);

      f.close();

    }
  }

  if(!topologyPath.size()){
    stringstream msg;
    msg 
      << "[Editor] No user topology specified!" << endl;
    dMsg(cerr, msg.str(), 1);
    exit(-1);
  }
  else{

    ifstream f(topologyPath.data(), ios::in);

    if(!f){
      stringstream msg;
      msg << "[Editor] Cannot read file '" << topologyPath << "'!" << endl;
      dMsg(cerr, msg.str(), 1);
      exit(-1);
    }
    else{
    
      string keyword;

      f >> keyword;
      if(keyword == "STF"){
        int     minNb = 0, maxNb = 0, vId = -1;

        f >> minNb;
        f >> maxNb;

        for(int i = 0; i < minNb + maxNb; i++){
          f >> keyword;
          f >> vId;

          if(keyword == "m")
            authorizedMinima_.push_back(vId);
          else if(keyword == "M")
            authorizedMaxima_.push_back(vId);
        }
      }

      f.close();
    }
  }

  initialMinValue_ = scalarField_.getMinimumValue();
  initialMaxValue_ = scalarField_.getMaximumValue();

  return 0;
}

int Editor::loadTriangleMesh(const string &path){

  triangleMesh_.setDebugLevel(debugLevel_);
  scalarField_.setDebugLevel(debugLevel_);

  ifstream f(path.data(), ios::in);

  if(!f){
    stringstream msg;
    msg << "[Editor] Cannot read file '" << path << "'." << endl;
    dMsg(cerr, msg.str(), 1);
    return -1;
  }

  triangleMesh_ << f;

  // rescale the mesh
  triangleMesh_.rescaleBoundingBox(9, 9, 9);

  scalarField_.registerTriangleMesh(&triangleMesh_);

  f.close();

  return 0;
}

int Editor::saveScalarField(const string &scalarPath, 
  const string &topologyPath){

  ofstream f(scalarPath.data(), ios::out);

  if(!f){
    stringstream msg;
    msg << "[Editor] Cannot write in file '" << scalarPath << "'." << endl;
    dMsg(cerr, msg.str(), 1);
    return -1;
  }

  scalarField_ >> f;

  f.close();

  // additionally save its topology
  scalarField_.saveTopology(topologyPath);

  stringstream msg;
  msg << "[Editor] Surface scalar field written in '" 
    << scalarPath << "'." << endl;
  msg << "[Editor] Field topology saved in '" << topologyPath << "'." << endl;
  dMsg(cout, msg.str(), 1);

  return 0;
}

