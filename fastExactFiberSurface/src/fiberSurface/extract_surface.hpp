/* 
 * file:                  vtkFiberSurfaceFilter.cpp
 * description:           extract fiber surface (inlined)
 * author(s):             Pavol Klacansky <pavol@klacansky.com>.
 * date:                  March 2015.
 */

#ifndef EXTRACT_SURFACE_HPP
#define EXTRACT_SURFACE_HPP

#include <cmath>

#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

/*
    marching tets cases
    n triangles + triangles (pairs of endpoints of edges to interpolate)
    all triangles counter clockwise
    inside black
*/
static const uint8_t mt[16][13] = {
    {}, // 0
    {1, 0, 1, 0, 2, 0, 3}, // 1
    {1, 0, 1, 1, 3, 1, 2}, // 2
    {2, 1, 2, 0, 2, 0, 3, 1, 2, 0, 3, 1, 3}, // 3
    {1, 0, 2, 1, 2, 2, 3}, // 4
    {2, 0, 1, 1, 2, 0, 3, 1, 2, 2, 3, 0, 3}, // 5
    {2, 0, 2, 0, 1, 1, 3, 0, 2, 1, 3, 2, 3}, // 6
    {1, 0, 3, 1, 3, 2, 3}, // 7
    {1, 0, 3, 2, 3, 1, 3}, // 8
    {2, 0, 1, 0, 2, 1, 3, 0, 2, 2, 3, 1, 3}, // 9
    {2, 1, 2, 0, 1, 0, 3, 1, 2, 0, 3, 2, 3}, // 10
    {1, 1, 2, 0, 2, 2, 3}, // 11
    {2, 0, 2, 1, 2, 0, 3, 1, 2, 1, 3, 0, 3}, // 12
    {1, 0, 1, 1, 2, 1, 3}, // 13
    {1, 0, 2, 0, 1, 0, 3}, // 14
    {} // 15
};

// stores for each triangle edge the adjacent face
// number of faces is determined from mt table
#define FACE_NULL UINT8_MAX
static const uint8_t mt_face[16][6] = {
    {},
    {3, 1, 2},
    {2, 0, 3},
    {3, 1, FACE_NULL, FACE_NULL, 2, 0},
    {3, 0, 1},
    {3, FACE_NULL, 2, 0, 1, FACE_NULL},
    {3, 2, FACE_NULL, FACE_NULL, 0, 1},
    {2, 0, 1},
    {1, 0, 2},
    {3, FACE_NULL, 2, 1, 0, FACE_NULL},
    {3, 2, FACE_NULL, FACE_NULL, 1, 0},
    {3, 1, 0},
    {3, FACE_NULL, 1, 0, 2, FACE_NULL},
    {3, 0, 2},
    {3, 2, 1},
    {}
};

/*
    left 1, right 2, middle 0
    0 grey, 1 white, 2 black
    first number denotes length
    if it is a pair of two vertices, first one is outside (e.g. clamped)
    0-2 vertices, 3-5 vertices as part of edges
*/
static const uint8_t mtr[27][18] = {
    {3, 0, 1, 2}, // carbon copy of old triangle, same as marching tetrahedra // 0
    {9, 3, 4, 2, 3, 5, 3, 4, 1, 2}, // 1
    {9, 3, 4, 2, 3, 5, 3, 4, 1, 2}, // 2
    {9, 0, 4, 3, 2, 4, 3, 4, 5, 2}, // 3
    {5, 3, 5, 4, 5, 2}, // 4
    {17, 3, 4, 3, 5, 4, 3, 4, 3, 3, 5, 4, 5, 3, 5, 2, 4, 5}, // 5
    {9, 0, 4, 3, 2, 4, 3, 4, 5, 2}, // 6
    {17, 3, 4, 4, 5, 3, 5, 3, 4, 4, 3, 4, 5, 3, 5, 4, 5, 2}, // 7
    {5, 3, 5, 4, 5, 2}, // 8
    {9, 0, 5, 4, 5, 3, 0, 1, 5, 4}, // 9
    {5, 3, 4, 1, 5, 4}, // 10
    {17, 3, 4, 5, 3, 3, 5, 3, 4, 1, 5, 4, 3, 4, 5, 4, 5, 3}, // 11
    {5, 0, 4, 3, 5, 3}, // 12
    {}, // 13 ALL NEGATIVE
    {12, 3, 4, 4, 3, 3, 5, 4, 3, 5, 3, 3, 5}, // 14
    {17, 0, 4, 3, 5, 3, 4, 3, 4, 5, 5, 4, 4, 3, 5, 4, 5, 3}, // 15
    {12, 3, 4, 4, 3, 4, 5, 3, 4, 4, 5, 5, 4}, // 16
    {12, 3, 5, 4, 5, 5, 3, 4, 5, 5, 4, 5, 3}, // 17
    {9, 0, 5, 4, 5, 3, 0, 1, 5, 4}, // 18
    {17, 3, 4, 5, 4, 3, 5, 3, 4, 1, 5, 4, 3, 5, 5, 4, 5, 3}, // 19
    {5, 3, 4, 1, 5, 4}, // 20
    {17, 0, 4, 3, 5, 3, 4, 3, 4, 5, 5, 3, 5, 3, 4, 5, 5, 4}, // 21
    {12, 3, 5, 4, 5, 5, 4, 3, 5, 5, 4, 5, 3}, // 22
    {12, 3, 4, 4, 3, 4, 5, 3, 4, 4, 5, 5, 4}, // 23
    {5, 0, 4, 3, 5, 3}, // 24
    {12, 3, 4, 4, 3, 5, 3, 3, 4, 5, 3, 3, 5}, // 25
    {} // 26 ALL POSITIVE
};


// should get inline so the stuck pressure is minimal
template<class dataType0, class dataType1> inline void extract_surface(
  const vtkIdType *tetIds,
  const float *pointSet,
  vtkPolyData *output, 
  const void *field0, 
  const void *field1,
  const double *o, 
  const double *n, 
  const double &distance, 
  const double &d_length, 
  const double &cur_length,
  const bool &visibleFibers,
  vtkFloatArray *textureCoordinates,
  vtkPoints *outputPoints,
  const int &edgeId,
  vtkIntArray *edgeIds){
  
  double h[4];
  float texCoords[2];
  vtkIdType triangle[3];

  // classify tet with respect to projected line segment
  dataType0 *fieldU = (dataType0 *) field0;
  dataType1 *fieldV = (dataType1 *) field1;
  uint8_t c = 0;
  
  for (int i = 0; i != 4; ++i) {
    h[i] = fieldU[tetIds[i]]*n[0] + fieldV[tetIds[i]]*n[1]
      - distance;
      c |= (h[i] <= 0.0f) << i;
      h[i] = std::abs(h[i]);
  }

  // process tet surface
  for (uint8_t i = 0; i != mt[c][0]; ++i) {
  
    double t[3];
    double tv[3][3];
    uint8_t tc = 0; // case of triangle

    for (int j = 0; j != 3; ++j) {
      const uint8_t *edge = mt[c] + 6*i + 2*j + 1;
      const double a = h[edge[0]]/(h[edge[0]] + h[edge[1]]);

      const float *p0 = &(pointSet[3*tetIds[edge[0]]]);
      const float *p1 = &(pointSet[3*tetIds[edge[1]]]);

      tv[j][0] = (1.0 - a)*p0[0] + a*p1[0];
      tv[j][1] = (1.0 - a)*p0[1] + a*p1[1];
      tv[j][2] = (1.0 - a)*p0[2] + a*p1[2];

      // parameter in range with respect to the projection
      const double v[] = {
        (1.0 - a)*fieldU[tetIds[edge[0]]] +
        a*fieldU[tetIds[edge[1]]] - o[0],
        (1.0 - a)*fieldV[tetIds[edge[0]]] +
        a*fieldV[tetIds[edge[1]]] - o[1],
      };
      t[j] = (v[0]*(-n[1]) + v[1]*n[0])/d_length;
      tc += ((t[j] < 0.0) + 2*(t[j] > 1.0))*((j == 2) ? 9 : 2*j + 1);
    }

    // subdivide the extracted tet surface
    for (uint8_t j = 0, counter = 0; j != mtr[tc][0]; ++j) {
      if (mtr[tc][j + 1] < 3) {
        triangle[counter++] = 
          outputPoints->InsertNextPoint(tv[mtr[tc][j + 1]]);
        
        if(visibleFibers){
          texCoords[0] = t[mtr[tc][j + 1]];
          texCoords[1] = texCoords[0]*d_length + cur_length;
        }
      } 
      else {
        ++j;
        const double a = 
          ((t[mtr[tc][j] - 3] < 0.0 ? 0.0 : 1.0) - t[mtr[tc][j] -3])/
            (t[mtr[tc][j + 1] - 3]  - t[mtr[tc][j] - 3]);
        triangle[counter++] = outputPoints->InsertNextPoint(
          (1.0 - a)*tv[mtr[tc][j] - 3][0] + a*tv[mtr[tc][j + 1] - 3][0],
          (1.0 - a)*tv[mtr[tc][j] - 3][1] + a*tv[mtr[tc][j + 1] - 3][1],
          (1.0 - a)*tv[mtr[tc][j] - 3][2] + a*tv[mtr[tc][j + 1] - 3][2]);
        if(visibleFibers){
          texCoords[0] = t[mtr[tc][j] - 3] < 0.0 ? 0.0f : 1.0f;
          texCoords[1] = texCoords[0]*d_length + cur_length;
        }
      }

      if(visibleFibers){
        textureCoordinates->InsertNextTuple(texCoords);
      }

      if (counter == 3) {
        if(edgeIds){
          edgeIds->InsertNextTuple1(edgeId);
        }
        output->InsertNextCell(VTK_TRIANGLE, 3, triangle);
        counter = 0;
      }
    }
  }
}

#endif
