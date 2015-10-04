/*
 * file:                main.cpp
 * description:         Example program that simplifies a scalar field given
 *                      some constraints on the output singularities.
 *                      Both input and output fields are stored as *wrl files
 *                      directly viewable with a standard VRML viewer or
 *                      Paraview, Blender, etc.
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

#include                  "Editor.h"

#include                  <fstream>
#include                  <iostream>
#include                  <sstream>
#include                  <string>

using namespace std;

int printUsage(char *binPath){

  cout << "Usage:" << endl;
  cout << binPath 
    << " -s <*off file> -f <scalar field file> -t <topology file> -r" << endl;

  cout << endl << "See README.txt for more information." << endl;

  return 0;
}

int main(int argc, char **argv){

  if(argc < 3){
    printUsage(argv[0]);
    return -1;
  }

  string  meshPath, scalarPath, topologyPath;

  for(int i = 0; i < argc; i++){

    if((string(argv[i]) == "-s")&&(i + 1 < argc)){
      meshPath = string(argv[i + 1]);
      i++;
    }
    if((string(argv[i]) == "-f")&&(i + 1 < argc)){
      scalarPath = string(argv[i + 1]);
      i++;
    }
    if((string(argv[i]) == "-t")&&(i + 1 < argc)){
      topologyPath = string(argv[i + 1]);
      i++;
    }
    if(string(argv[i]) == "-h"){
      printUsage(argv[0]);
      return -1;
    }
  }

  if(!meshPath.size()){
    printUsage(argv[0]);
    return -1;
  }

  Editor editor;

  editor.loadTriangleMesh(meshPath);
  editor.loadScalarField(scalarPath, topologyPath);
  editor.exportToVRML("input.wrl");
 
  editor.simplifyField();

  editor.saveScalarField("outputField.ssf", "outputTopology.sft");
  editor.exportToVRML("output.wrl");

  return 0;
}
