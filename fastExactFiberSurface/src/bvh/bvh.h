/* 
 * file:                  bvh.h
 * description:           Bounding Volume Hierarchy.
 * author(s):             Pavol Klacansky <pavol@klacansky.com>.
 * date:                  March 2015.
 */

/*
	PURPOSE: Acceleration data structure for testing segment with aabb

	REFERENCE: Shallow Bounding Volume Hierarchies for Fast SIMD Ray Tracing of
	           Incoherent Rays
*/

#ifndef BVH_H
#define BVH_H

#include <stdint.h>

// TODO add support for width 4 (SSE) and 8 (AVX) - probably in separate files
// MUST BE POWER OF TWO
#define BVH_WIDTH 8

// bbox [min_x,min_y,...,max_x,max_y,...]
struct BvhNode {
	float bbox[4][BVH_WIDTH]; // 16*BVH_WIDTH bytes
	int child[BVH_WIDTH]; // 4*BVH_WIDTH bytes, INT_MIN = empty leaf, negative = leaf
};

struct Bvh {
	struct BvhNode *nodes;
	uint32_t length;
	uint32_t capacity;

	uint32_t *indices;
	uint32_t indices_length;
	uint8_t min_cells_per_leaf;
};

struct Bvh *
bvh_alloc(const uint32_t size,
          const uint8_t min_cells_per_leaf);

void
bvh_build(struct Bvh *bvh, const float *bboxs_orig);

void
bvh_free(struct Bvh *bvh);

uint32_t
bvh_intersect(const struct Bvh *bvh,
              uint32_t *buffer,
              const uint32_t size,
              const float *line);

void
bvh_print(const struct Bvh *bvh);

void
bvh_stats(const struct Bvh *bvh);

#endif
