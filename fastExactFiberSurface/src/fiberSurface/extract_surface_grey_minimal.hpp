/*
 * file:                  vtkFiberSurfaceFilter.cpp
 * description:           extract minimal fiber surface with grey cases (inlined)
 * author(s):             Pavol Klacansky <pavol@klacansky.com>.
 * date:                  March 2015.
 */

#ifndef EXTRACT_SURFACE_GREY_MINIMAL_HPP
#define EXTRACT_SURFACE_GREY_MINIMAL_HPP

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
static const uint8_t fs_mt[][9] = {
    {},
    {3, 1, 3, 2},
    {3, 1, 2, 3},
    {3, 3, 0, 2},
    {},
    {4, 2, 3, 4, 5},
    {3, 3, 2, 0},
    {4, 3, 2, 5, 4},
    {},
    {3, 3, 1, 0},
    {},
    {4, 3, 1, 4, 6},
    {},
    {},
    {5, 4, 6, 3, 4, 5},
    {4, 0, 3, 5, 6},
    {5, 5, 4, 3, 5, 6},
    {5, 3, 5, 6, 4, 6},
    {3, 1, 3, 0},
    {4, 1, 3, 6, 4},
    {},
    {4, 3, 0, 6, 5},
    {5, 6, 5, 3, 6, 4},
    {5, 3, 4, 5, 6, 5},
    {},
    {5, 3, 6, 4, 5, 4},
    {},
    {3, 0, 1, 2},
    {},
    {4, 1, 2, 4, 7},
    {},
    {},
    {5, 4, 5, 2, 4, 7},
    {4, 2, 0, 5, 7},
    {5, 5, 7, 2, 5, 4},
    {5, 2, 4, 7, 5, 7},
    {},
    {},
    {5, 4, 7, 1, 4, 6},
    {},
    {},
    {6, 4, 5, 4, 6, 4, 7},
    {5, 5, 6, 0, 5, 7},
    {6, 5, 7, 5, 6, 5, 4},
    {8, 5, 6, 4, 6, 4, 7, 5, 7},
    {4, 0, 1, 6, 7},
    {5, 6, 4, 1, 6, 7},
    {5, 1, 6, 7, 4, 7},
    {5, 6, 7, 0, 6, 5},
    {6, 6, 5, 6, 7, 6, 4},
    {8, 6, 7, 4, 7, 4, 5, 6, 5},
    {5, 0, 5, 7, 6, 7},
    {8, 6, 4, 5, 4, 5, 7, 6, 7},
    {6, 4, 7, 5, 7, 6, 7},
    {3, 0, 2, 1},
    {4, 2, 1, 7, 4},
    {},
    {4, 0, 2, 7, 5},
    {5, 7, 4, 2, 7, 5},
    {5, 2, 7, 5, 4, 5},
    {},
    {5, 2, 5, 4, 7, 4},
    {},
    {4, 1, 0, 7, 6},
    {5, 7, 6, 1, 7, 4},
    {5, 1, 4, 6, 7, 6},
    {5, 7, 5, 0, 7, 6},
    {6, 7, 4, 7, 6, 7, 5},
    {8, 4, 6, 7, 6, 7, 5, 4, 5},
    {5, 0, 7, 6, 5, 6},
    {8, 7, 6, 5, 6, 5, 4, 7, 4},
    {6, 7, 6, 5, 6, 4, 6},
    {},
    {5, 1, 7, 4, 6, 4},
    {},
    {5, 0, 6, 5, 7, 5},
    {8, 6, 5, 7, 5, 7, 4, 6, 4},
    {6, 7, 5, 4, 5, 6, 5},
    {},
    {6, 5, 4, 7, 4, 6, 4},
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
const uint8_t fs_mtr[][18] = {
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

// 81 cases for a fiber surface rectangle
const uint8_t fs_mre[][23] = {
    {6, 3, 0, 1, 3, 1, 2},
    {13, 1, 4, 7, 4, 5, 1, 2, 4, 7, 2, 3, 4, 7},
    {12, 1, 2, 3, 1, 3, 4, 7, 1, 4, 7, 4, 5},
    {13, 2, 5, 4, 5, 6, 2, 3, 5, 4, 3, 0, 5, 4},
    {9, 2, 3, 5, 6, 3, 4, 7, 5, 6},
    {20, 2, 5, 4, 5, 6, 2, 3, 5, 4, 3, 4, 7, 5, 4, 4, 7, 4, 5, 5, 4},
    {12, 2, 3, 0, 2, 0, 5, 4, 2, 5, 4, 5, 6},
    {20, 5, 6, 4, 5, 5, 4, 2, 4, 5, 5, 6, 2, 3, 4, 5, 3, 4, 7, 4, 5},
    {9, 2, 3, 4, 7, 2, 4, 7, 5, 6},
    {13, 3, 6, 5, 6, 7, 3, 0, 6, 5, 0, 1, 6, 5},
    {},
    {22, 3, 4, 7, 6, 7, 4, 7, 4, 5, 6, 5, 4, 7, 6, 5, 6, 7, 4, 5, 1, 6, 5},
    {9, 3, 0, 6, 7, 0, 5, 4, 6, 7},
    {5, 3, 4, 7, 6, 7},
    {17, 3, 4, 7, 6, 7, 4, 7, 4, 5, 5, 4, 4, 7, 5, 4, 6, 7},
    {20, 3, 6, 5, 6, 7, 3, 0, 6, 5, 0, 5, 4, 6, 5, 5, 4, 5, 6, 6, 5},
    {},
    {17, 3, 4, 7, 6, 7, 4, 7, 5, 6, 6, 7, 5, 6, 6, 5, 6, 7},
    {12, 3, 0, 1, 3, 1, 6, 5, 3, 6, 5, 6, 7},
    {22, 1, 6, 5, 4, 5, 6, 5, 6, 7, 4, 7, 6, 5, 4, 7, 4, 5, 6, 7, 3, 4, 7},
    {},
    {20, 6, 7, 5, 6, 6, 5, 3, 5, 6, 6, 7, 3, 0, 5, 6, 0, 5, 4, 5, 6},
    {17, 6, 7, 5, 6, 6, 5, 6, 7, 4, 7, 5, 6, 6, 7, 3, 4, 7},
    {},
    {9, 3, 0, 5, 4, 3, 5, 4, 6, 7},
    {17, 5, 4, 6, 7, 4, 7, 5, 4, 4, 7, 4, 5, 6, 7, 3, 4, 7},
    {5, 4, 7, 6, 7, 3},
    {13, 0, 7, 6, 7, 4, 0, 1, 7, 6, 1, 2, 7, 6},
    {9, 1, 2, 4, 5, 2, 7, 6, 4, 5},
    {20, 4, 5, 7, 4, 4, 7, 1, 7, 4, 4, 5, 1, 2, 7, 4, 2, 7, 6, 7, 4},
    {},
    {5, 2, 7, 6, 5, 6},
    {},
    {22, 0, 5, 4, 7, 4, 5, 4, 5, 6, 7, 6, 5, 4, 7, 6, 7, 4, 5, 6, 2, 7, 6},
    {17, 5, 6, 4, 5, 5, 4, 5, 6, 7, 6, 4, 5, 5, 6, 2, 7, 6},
    {17, 4, 7, 5, 6, 7, 6, 4, 7, 7, 6, 7, 4, 5, 6, 2, 7, 6},
    {9, 0, 1, 7, 4, 1, 6, 5, 7, 4},
    {5, 1, 6, 5, 4, 5},
    {17, 4, 5, 7, 4, 4, 7, 4, 5, 6, 5, 7, 4, 4, 5, 1, 6, 5},
    {5, 0, 5, 4, 7, 4},
    {},
    {12, 4, 5, 7, 4, 4, 7, 4, 5, 5, 4, 7, 4},
    {17, 0, 5, 4, 7, 4, 5, 4, 5, 6, 6, 5, 5, 4, 6, 5, 7, 4},
    {12, 5, 6, 4, 5, 5, 4, 5, 6, 6, 5, 4, 5},
    {12, 4, 7, 5, 6, 7, 4, 5, 6, 6, 5, 7, 4},
    {20, 0, 7, 6, 7, 4, 0, 1, 7, 6, 1, 6, 5, 7, 6, 6, 5, 6, 7, 7, 6},
    {17, 1, 6, 5, 4, 5, 6, 5, 6, 7, 7, 6, 6, 5, 7, 6, 4, 5},
    {},
    {},
    {12, 6, 7, 5, 6, 6, 5, 6, 7, 7, 6, 5, 6},
    {},
    {17, 0, 5, 4, 7, 4, 5, 4, 6, 7, 7, 4, 6, 7, 7, 6, 7, 4},
    {12, 5, 4, 6, 7, 4, 5, 6, 7, 7, 6, 4, 5},
    {12, 4, 7, 6, 7, 7, 4, 7, 4, 6, 7, 7, 6},
    {12, 0, 1, 2, 0, 2, 7, 6, 0, 7, 6, 7, 4},
    {20, 1, 4, 7, 4, 5, 1, 2, 4, 7, 2, 7, 6, 4, 7, 7, 6, 7, 4, 4, 7},
    {9, 1, 2, 7, 6, 1, 7, 6, 4, 5},
    {22, 2, 7, 6, 5, 6, 7, 6, 7, 4, 5, 4, 7, 6, 5, 4, 5, 6, 7, 4, 0, 5, 4},
    {17, 2, 7, 6, 5, 6, 7, 6, 7, 4, 4, 7, 7, 6, 4, 7, 5, 6},
    {17, 2, 7, 6, 5, 6, 7, 6, 4, 5, 5, 6, 4, 5, 5, 4, 5, 6},
    {},
    {},
    {5, 7, 6, 5, 6, 2},
    {20, 7, 4, 6, 7, 7, 6, 0, 6, 7, 7, 4, 0, 1, 6, 7, 1, 6, 5, 6, 7},
    {},
    {17, 7, 6, 4, 5, 6, 5, 7, 6, 6, 5, 6, 7, 4, 5, 1, 6, 5},
    {17, 7, 4, 6, 7, 7, 6, 7, 4, 5, 4, 6, 7, 7, 4, 0, 5, 4},
    {12, 7, 4, 6, 7, 7, 6, 7, 4, 4, 7, 6, 7},
    {12, 7, 6, 4, 5, 6, 7, 4, 5, 5, 4, 6, 7},
    {},
    {},
    {12, 7, 6, 5, 6, 6, 7, 6, 7, 5, 6, 6, 5},
    {9, 0, 1, 6, 5, 0, 6, 5, 7, 4},
    {17, 1, 6, 5, 4, 5, 6, 5, 7, 4, 4, 5, 7, 4, 4, 7, 4, 5},
    {5, 6, 5, 4, 5, 1},
    {17, 6, 5, 7, 4, 5, 4, 6, 5, 5, 4, 5, 6, 7, 4, 0, 5, 4},
    {12, 6, 5, 7, 4, 5, 6, 7, 4, 4, 7, 5, 6},
    {12, 6, 5, 4, 5, 5, 6, 5, 6, 4, 5, 5, 4},
    {5, 5, 4, 7, 4, 0},
    {12, 5, 4, 7, 4, 4, 5, 4, 5, 7, 4, 4, 7},
    {}
};



template<typename dataType0, typename dataType1> inline void
extract_surface_grey_minimal(const vtkIdType *tetIds,
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
	double p[4][3];
	double t[4];
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
	if (!fs_mt[c][0])
		return;

	// process tet surface
	uint8_t tc = 0; // case for triangle or quad
	for (uint8_t i = 0, of = 0, shift = 1; i != fs_mt[c][0]; ++i, ++of, shift *= 3) {
		const uint8_t *edge = fs_mt[c] + i + 1;
		double ft[2];

		// grey vertex
		if (edge[0] < 4) {
			p[of][0] = pointSet[3*tetIds[edge[0]]];
			p[of][1] = pointSet[3*tetIds[edge[0]] + 1];
			p[of][2] = pointSet[3*tetIds[edge[0]] + 2];

			ft[0] = fieldU[tetIds[edge[0]]];
			ft[1] = fieldV[tetIds[edge[0]]];

		// edge
		} else {
			++i; // reading two values

			// apparently C++ does not like narrowing
			uint8_t e[2];
			e[0] = edge[0] - 4;
			e[1] = edge[1] - 4;
			const double a = h[e[0]]/(h[e[0]] + h[e[1]]);

			const float *p0 = pointSet + 3*tetIds[e[0]];
			const float *p1 = pointSet + 3*tetIds[e[1]];
			p[of][0] = (1.0 - a)*p0[0] + a*p1[0];
			p[of][1] = (1.0 - a)*p0[1] + a*p1[1];
			p[of][2] = (1.0 - a)*p0[2] + a*p1[2];

			ft[0] = (1.0 - a)*fieldU[tetIds[e[0]]] + a*fieldU[tetIds[e[1]]];
			ft[1] = (1.0 - a)*fieldV[tetIds[e[0]]] + a*fieldV[tetIds[e[1]]];
		}
		// project onto a polyline segment
		t[of] = ((ft[0] - o[0])*-n[1] + (ft[1] - o[1])*n[0])*length_inv;

		/*
			classify
			t < 0 -> minus (value 1)
			t >= 0 &&  t <= 1 -> grey (value 0)
			t > 1 -> plus (value 2)
		*/
		tc += ((t[of] < 0.0) + 2*(t[of] > 1.0))*shift;
	}

	// decide if we need to process quad or triangle
	const uint8_t *table;
	uint8_t n_vertices;
	if (fs_mt[c][0] > 6) {
		table = fs_mre[tc];
		n_vertices = 4;
	} else {
		table = fs_mtr[tc];
		n_vertices = 3;
	}

	// process surface to obtain fiber surface
	for (uint8_t i = 0, counter = 0; i != table[0]; ++i) {
		// edge (n_vertices + v0, ...)
		if (table[i + 1] >= n_vertices) {
			++i; // reading 2 values

			// trust me C++
			uint8_t offs[2];
			offs[0] = table[i] - n_vertices;
			offs[1] = table[i + 1] - n_vertices;

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
			triangle[counter++] = outputPoints->InsertNextPoint(p[table[i + 1]]);

			if (visibleFibers) {
				texCoords[0] = t[table[i + 1]];
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

#endif
