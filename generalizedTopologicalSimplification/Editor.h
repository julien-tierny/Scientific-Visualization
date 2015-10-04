/*
 * file:                Editor.h
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

#ifndef                 _EDITOR_H
#define                 _EDITOR_H

#include                <SurfaceScalarField.h>

using namespace std;

class Editor : public Debug{

  public:

    Editor(){};
    ~Editor(){};

    int exportToVRML(const string &path);

    int loadScalarField(const string &scalarPath, const string &topologyPath);

    int loadTriangleMesh(const string &path);
  
    int saveScalarField(const string &scalarPath, const string &topologyPath);

    inline int simplifyField() { 
      return scalarField_.simplify(authorizedMinima_, authorizedMaxima_);};

  protected:

    vector<int>             authorizedMinima_, authorizedMaxima_;
    double                  initialMinValue_, initialMaxValue_;
    SurfaceScalarField      scalarField_;
    TriangleMesh            triangleMesh_;
};
#endif
