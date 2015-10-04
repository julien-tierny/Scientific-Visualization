/*
 * file:                  vtkFiberSurfaceFilter.cpp
 * description:           extract minimal fiber surface with grey cases (inlined)
 * author(s):             Pavol Klacansky <pavol@klacansky.com>.
 * date:                  March 2015.
 */

#ifndef EXTRACT_SURFACE_GREY_HPP
#define EXTRACT_SURFACE_GREY_HPP

#include <cmath>

#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

/*
    0 = grey (on surface)
    1 = white (outside)
    2 = black (inside)
*/

/*
    tet indexing

    3 --- 5 --- 2 --- 5 --- 3
     \         / \         /
      \       /   \       /
       3  1  1  3  0  0  4
        \   /       \   /
         \ /         \ /
          0 --- 2 --- 1
           \         /
            \       /
             3  2  4
              \   /
               \ /
                3
*/
// marching tetrahedra surface with grey vertices
// TODO grey vertex indices are listed in pairs, this is because of the original
//      implementation using SIMD to avoid branching; edge pair has distinct values (e.g. 0, 2)
static const uint8_t g_mt[81][13] = {
    {},
    {6, 1, 1, 3, 3, 2, 2},
    {6, 1, 1, 2, 2, 3, 3},
    {6, 3, 3, 0, 0, 2, 2},
    {},
    {6, 2, 2, 3, 3, 0, 1},
    {6, 3, 3, 2, 2, 0, 0},
    {6, 3, 3, 2, 2, 1, 0},
    {},
    {6, 3, 3, 1, 1, 0, 0},
    {},
    {6, 3, 3, 1, 1, 0, 2},
    {},
    {},
    {6, 0, 2, 3, 3, 0, 1},
    {6, 0, 0, 3, 3, 1, 2},
    {6, 1, 0, 3, 3, 1, 2},
    {6, 3, 3, 1, 2, 0, 2},
    {6, 1, 1, 3, 3, 0, 0},
    {6, 1, 1, 3, 3, 2, 0},
    {},
    {6, 3, 3, 0, 0, 2, 1},
    {6, 2, 1, 3, 3, 2, 0},
    {6, 3, 3, 0, 1, 2, 1},
    {},
    {6, 3, 3, 2, 0, 1, 0},
    {},
    {6, 0, 0, 1, 1, 2, 2},
    {},
    {6, 1, 1, 2, 2, 0, 3},
    {},
    {},
    {6, 0, 1, 2, 2, 0, 3},
    {6, 2, 2, 0, 0, 1, 3},
    {6, 1, 3, 2, 2, 1, 0},
    {6, 2, 2, 0, 3, 1, 3},
    {},
    {},
    {6, 0, 3, 1, 1, 0, 2},
    {},
    {},
    {6, 0, 1, 0, 2, 0, 3},
    {6, 1, 2, 0, 0, 1, 3},
    {6, 1, 3, 1, 2, 1, 0},
    {12, 1, 2, 0, 2, 0, 3, 1, 2, 0, 3, 1, 3},
    {6, 0, 0, 1, 1, 2, 3},
    {6, 2, 0, 1, 1, 2, 3},
    {6, 1, 1, 2, 3, 0, 3},
    {6, 2, 3, 0, 0, 2, 1},
    {6, 2, 1, 2, 3, 2, 0},
    {12, 2, 3, 0, 3, 0, 1, 2, 3, 0, 1, 2, 1},
    {6, 0, 0, 1, 3, 2, 3},
    {12, 2, 0, 1, 0, 1, 3, 2, 0, 1, 3, 2, 3},
    {6, 0, 3, 1, 3, 2, 3},
    {6, 0, 0, 2, 2, 1, 1},
    {6, 2, 2, 1, 1, 3, 0},
    {},
    {6, 0, 0, 2, 2, 3, 1},
    {6, 3, 0, 2, 2, 3, 1},
    {6, 2, 2, 3, 1, 0, 1},
    {},
    {6, 2, 2, 1, 0, 3, 0},
    {},
    {6, 1, 1, 0, 0, 3, 2},
    {6, 3, 2, 1, 1, 3, 0},
    {6, 1, 1, 0, 2, 3, 2},
    {6, 3, 1, 0, 0, 3, 2},
    {6, 3, 0, 3, 2, 3, 1},
    {12, 0, 2, 3, 2, 3, 1, 0, 2, 3, 1, 0, 1},
    {6, 0, 0, 3, 2, 1, 2},
    {12, 3, 2, 1, 2, 1, 0, 3, 2, 1, 0, 3, 0},
    {6, 3, 2, 1, 2, 0, 2},
    {},
    {6, 1, 1, 3, 0, 2, 0},
    {},
    {6, 0, 0, 2, 1, 3, 1},
    {12, 2, 1, 3, 1, 3, 0, 2, 1, 3, 0, 2, 0},
    {6, 3, 1, 0, 1, 2, 1},
    {},
    {6, 1, 0, 3, 0, 2, 0},
    {}
};

/*
    triangle numbering

        2
       / \
      1   0
     /     \
    0 - 2 - 1
*/
// 27 cases for a fiber surface triangle
const uint8_t g_mtr[27][18] = {
    {3, 2, 0, 1},
    {9, 1, 2, 3, 5, 1, 3, 5, 3, 4},
    {9, 1, 2, 3, 5, 1, 3, 5, 3, 4},
    {9, 2, 0, 4, 3, 2, 4, 3, 4, 5},
    {5, 2, 3, 5, 4, 5},
    {17, 2, 3, 5, 4, 5, 3, 5, 3, 4, 4, 5, 3, 4, 4, 3, 4, 5},
    {9, 2, 0, 4, 3, 2, 4, 3, 4, 5},
    {17, 4, 5, 3, 4, 4, 3, 4, 5, 3, 5, 3, 4, 4, 5, 2, 3, 5},
    {5, 3, 5, 4, 5, 2},
    {9, 0, 1, 5, 4, 0, 5, 4, 5, 3},
    {5, 1, 5, 4, 3, 4},
    {17, 3, 4, 5, 3, 3, 5, 3, 4, 5, 4, 5, 3, 3, 4, 1, 5, 4},
    {5, 0, 4, 3, 5, 3},
    {},
    {12, 3, 4, 4, 3, 5, 3, 3, 4, 5, 3, 3, 5},
    {17, 0, 4, 3, 5, 3, 4, 3, 4, 5, 5, 3, 4, 5, 5, 4, 5, 3},
    {12, 4, 5, 5, 4, 3, 4, 4, 5, 3, 4, 4, 3},
    {12, 3, 5, 4, 5, 5, 4, 3, 5, 5, 4, 5, 3},
    {9, 0, 1, 5, 4, 0, 5, 4, 5, 3},
    {17, 1, 5, 4, 3, 4, 5, 4, 5, 3, 3, 4, 5, 3, 3, 5, 3, 4},
    {5, 5, 4, 3, 4, 1},
    {17, 5, 3, 4, 5, 5, 4, 5, 3, 4, 3, 4, 5, 5, 3, 0, 4, 3},
    {12, 5, 3, 3, 5, 4, 5, 5, 3, 4, 5, 5, 4},
    {12, 5, 4, 3, 4, 4, 3, 5, 4, 4, 3, 4, 5},
    {5, 4, 3, 5, 3, 0},
    {12, 4, 3, 5, 3, 3, 5, 4, 3, 3, 5, 3, 4},
    {}
};


template<typename dataType0, typename dataType1> inline void
extract_surface_grey(const vtkIdType *tetIds,
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
                     vtkIntArray *edgeIds)
{
	const double length_inv = 1.0/d_length;
	double h[4];
	double p[3][3];
	double t[3];
	float texCoords[2];
	vtkIdType triangle[3];

	// classify tet with respect to projected line segment
	const dataType0 *fieldU = static_cast<const dataType0 *>(field0);
	const dataType1 *fieldV = static_cast<const dataType1 *>(field1);
	uint8_t c = 0;

	// if we assumed floats, we could precompute derived field and save bandwidth
	for (int i = 0, shift = 1; i != 4; ++i, shift *= 3) {
		h[i] = fieldU[tetIds[i]]*n[0] + fieldV[tetIds[i]]*n[1] - distance;
		c += ((h[i] < 0.0) + 2*(h[i] > 0.0))*shift;
		h[i] = std::abs(h[i]);
	}

	// no surface, terminate
	if (!g_mt[c][0])
		return;

	// process tet surface
	for (uint8_t i = 0; i != g_mt[c][0]; i += 6) {
		uint8_t tc = 0; // case for triangle

		// process triangle
		for (int j = 0, shift = 1; j != 3; ++j, shift *= 3) {
			const uint8_t *edge = g_mt[c] + i + 2*j + 1;
			double ft[2];

			// grey vertex
			if (edge[0] == edge[1]) {
				p[j][0] = pointSet[3*tetIds[edge[0]]];
				p[j][1] = pointSet[3*tetIds[edge[0]] + 1];
				p[j][2] = pointSet[3*tetIds[edge[0]] + 2];

				ft[0] = fieldU[tetIds[edge[0]]];
				ft[1] = fieldV[tetIds[edge[0]]];

			// edge
			} else {
				// apparently C++ does not like narrowing
				const double a = h[edge[0]]/(h[edge[0]] + h[edge[1]]);

				const float *p0 = pointSet + 3*tetIds[edge[0]];
				const float *p1 = pointSet + 3*tetIds[edge[1]];
				p[j][0] = (1.0 - a)*p0[0] + a*p1[0];
				p[j][1] = (1.0 - a)*p0[1] + a*p1[1];
				p[j][2] = (1.0 - a)*p0[2] + a*p1[2];

				ft[0] = (1.0 - a)*fieldU[tetIds[edge[0]]] + a*fieldU[tetIds[edge[1]]];
				ft[1] = (1.0 - a)*fieldV[tetIds[edge[0]]] + a*fieldV[tetIds[edge[1]]];
			}
			// project onto a polyline segment
			t[j] = ((ft[0] - o[0])*-n[1] + (ft[1] - o[1])*n[0])*length_inv;

			/*
				classify
				t < 0 -> minus (value 1)
				t >= 0 &&  t <= 1 -> grey (value 0)
				t > 1 -> plus (value 2)
			*/
			tc += ((t[j] < 0.0) + 2*(t[j] > 1.0))*shift;
		}

		// process surface to obtain fiber surface
		for (uint8_t j = 0, counter = 0; j != g_mtr[tc][0]; ++j) {
			// edge (n_vertices + v0, ...)
			if (g_mtr[tc][j + 1] >= 3) {
				++j; // reading 2 values

				// trust me C++
				uint8_t offs[2];
				offs[0] = g_mtr[tc][j] - 3;
				offs[1] = g_mtr[tc][j + 1] - 3;

				// clamp the vertex either to 0 or 1
				const double a = ((t[offs[0]] < 0.0 ? 0.0 : 1.0) - t[offs[0]])/
								 (t[offs[1]] - t[offs[0]]);

				triangle[counter++] = outputPoints->InsertNextPoint(
					(1.0 - a)*p[offs[0]][0] + a*p[offs[1]][0],
					(1.0 - a)*p[offs[0]][1] + a*p[offs[1]][1],
					(1.0 - a)*p[offs[0]][2] + a*p[offs[1]][2]
				);

				if (visibleFibers) {
					texCoords[0] = t[offs[0]] < 0.0 ? 0.0 : 1.0;
					texCoords[1] = texCoords[0]*d_length + cur_length;
				}

			// vertex
			} else {
				triangle[counter++] = outputPoints->InsertNextPoint(p[g_mtr[tc][j + 1]]);

				if (visibleFibers) {
					texCoords[0] = t[g_mtr[tc][j + 1]];
					texCoords[1] = texCoords[0]*d_length + cur_length;
				}
			}

			if (visibleFibers)
				textureCoordinates->InsertNextTuple(texCoords);

			if (counter == 3) {
				if(edgeIds)
					edgeIds->InsertNextTuple1(edgeId);

				output->InsertNextCell(VTK_TRIANGLE, 3, triangle);
				counter = 0;
			}
		}
	}
}

#endif
