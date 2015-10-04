/* 
 * file:                  bvh.c
 * description:           Bounding Volume Hierarchy.
 * author(s):             Pavol Klacansky <pavol@klacansky.com>.
 * date:                  March 2015.
 */

#include <assert.h>
#include <float.h>
#include <inttypes.h>
#include <immintrin.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time_.h>


#include "bvh.h"


#define STACK_LENGTH 128

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))


static const int leaf_empty = INT_MIN;
static const uint32_t offset_mask = ~0U >> 5;

enum Axis {AXIS_X, AXIS_Y, N_AXES};


// 31-27 bits (4) represents number of children, starting from 1 to 16
static uint8_t
child_length(const uint32_t child)
{
	return (child >> 27 & 0x0F) + 1;
}


// TODO alignment of node bboxes
// register pressure should be reduced by compiler, we prefer readability
static void
intersect_8(float *restrict isect, const float *restrict bboxs, const float *restrict o, const float *restrict d_inv)
{
	__m256 ymm0 = _mm256_loadu_ps(bboxs);
	__m256 ymm1 = _mm256_loadu_ps(bboxs + 2*BVH_WIDTH);
	// subtract origin for first dim
	__m256 ymm2 = _mm256_broadcast_ss(o);
	ymm0 = _mm256_sub_ps(ymm0, ymm2);
	ymm1 = _mm256_sub_ps(ymm1, ymm2);

	// divide by direction for first dim
	ymm2 = _mm256_broadcast_ss(d_inv);
	ymm0 = _mm256_mul_ps(ymm0, ymm2);
	ymm1 = _mm256_mul_ps(ymm1, ymm2);

	// compute min/max for first dim
	__m256 ymm3 = _mm256_min_ps(ymm0, ymm1);
	__m256 ymm4 = _mm256_max_ps(ymm0, ymm1);


	// process second dimension
	ymm0 = _mm256_loadu_ps(bboxs + BVH_WIDTH);
	ymm1 = _mm256_loadu_ps(bboxs + 3*BVH_WIDTH);

	// miny - oy, maxy - oy
	ymm2 = _mm256_broadcast_ss(o + 1);
	ymm0 = _mm256_sub_ps(ymm0, ymm2);
	ymm1 = _mm256_sub_ps(ymm1, ymm2);

	// divide by dy
	ymm2 = _mm256_broadcast_ss(d_inv + 1);
	ymm0 = _mm256_mul_ps(ymm0, ymm2);
	ymm1 = _mm256_mul_ps(ymm1, ymm2);

	// min/max for second dim
	__m256 ymm5 = _mm256_min_ps(ymm0, ymm1);
	__m256 ymm6 = _mm256_max_ps(ymm0, ymm1);

	// total min/max
	ymm0 = _mm256_max_ps(ymm3, ymm5);
	ymm1 = _mm256_min_ps(ymm4, ymm6);

	// compute if the box is intersected
	ymm2 = _mm256_cmp_ps(ymm0, ymm1, _CMP_LT_OQ);
	ymm3 = _mm256_cmp_ps(ymm0, _mm256_set1_ps(1.0f), _CMP_LT_OQ);
	ymm2 = _mm256_and_ps(ymm2, ymm3);
	ymm3 = _mm256_cmp_ps(ymm1, _mm256_set1_ps(0.0f), _CMP_GT_OQ);
	ymm2 = _mm256_and_ps(ymm2, ymm3);

	_mm256_store_ps(isect, ymm2);
}

static void
intersect(int *restrict isect, const float *restrict bbox, const float *restrict o, const float *restrict d_inv, const int i)
{
	float tmin[2], tmax[2];

	// origin inside
	tmin[0] = (bbox[i] - o[0])*d_inv[0];
	tmax[0] = (bbox[2*BVH_WIDTH + i] - o[0])*d_inv[0];
	if (tmin[0] > tmax[0]) {
		const float tmp = tmin[0];
		tmin[0] = tmax[0];
		tmax[0] = tmp;
	}

	tmin[1] = (bbox[BVH_WIDTH + i] - o[1])*d_inv[1];
	tmax[1] = (bbox[3*BVH_WIDTH + i] - o[1])*d_inv[1];
	if (tmin[1] > tmax[1]) {
        const float tmp = tmin[1];
        tmin[1] = tmax[1];
        tmax[1] = tmp;
    }

	tmin[0] = (tmin[0] > tmin[1]) ? tmin[0] : tmin[1];
	tmax[0] = (tmax[0] < tmax[1]) ? tmax[0] : tmax[1];

	isect[i] = tmin[0] <= tmax[0] && tmin[0] < 1.0f && tmax[0] > 0.0f;
}


static void
bbox_swap(float *a, float *b)
{
	__m128 xmm0 = _mm_load_ps(a);
	__m128 xmm1 = _mm_load_ps(b);
	_mm_store_ps(a, xmm1);
	_mm_store_ps(b, xmm0);
}


static void
bbox_update(float *dest, const float *src)
{
	dest[0] = MIN(dest[0], src[0]);
	dest[1] = MIN(dest[1], src[1]);
	dest[2] = MAX(dest[2], src[2]);
    dest[3] = MAX(dest[3], src[3]);
}


static uint32_t
split(float *restrict bbox, // two new bounding boxes
      float *restrict bboxs, // input array of bounding boxes to split
      uint32_t *restrict indices, // index array with cell ids
      const uint32_t length, // number of bounding boxes in the array
      const enum Axis axis, // axis along which to split
      const float mid) // value along which to split
{
	bbox[0] = bbox[1] = bbox[4] = bbox[5] = FLT_MAX;
	bbox[2] = bbox[3] = bbox[6] = bbox[7] = -FLT_MAX;
	uint32_t offset = 0;

	for (uint32_t i = 0; i != length; ++i) {
		const float *cur = bboxs + 4*i;
		// TODO SIMD
		// max of the axis is smaller than mid
		if (cur[2 + axis] <= mid) {
			// most go first because pointer aliasing
			bbox_update(bbox, cur);

			// swap
			bbox_swap(bboxs + 4*offset, bboxs + 4*i);

			const uint32_t tmp = indices[offset];
			indices[offset] = indices[i];
			indices[i] = tmp;

			++offset;
		} else
			bbox_update(bbox + 4, cur);
	}

	// balance the tree
	// nothing on the min side (this happens only for smallish buckets)
	if (offset == 0 || offset == length) {
		offset = length/2;
		// recompute bboxes
		bbox[0] = bbox[1] = bbox[4] = bbox[5] = FLT_MAX;
	    bbox[2] = bbox[3] = bbox[6] = bbox[7] = -FLT_MAX;
		for (uint32_t i = 0; i != offset; ++i)
			bbox_update(bbox, bboxs + 4*i);
		const uint32_t n = length - offset;
		for (uint32_t i = 0; i != n; ++i)
			bbox_update(bbox + 4, bboxs + 4*(offset + i));
	}

	return offset;
}


struct Bvh *
bvh_alloc(const uint32_t size, const uint8_t min_cells_per_leaf)
{
	if (min_cells_per_leaf > 16) {
		fprintf(stderr, "BVH: Maximum number of cells per leaf is 16\n");
		return NULL;
	}

	struct Bvh *qbvh = malloc(sizeof(struct Bvh));
	qbvh->nodes = malloc(1024*sizeof(struct BvhNode));
	qbvh->length = 0;
	qbvh->capacity = 1024;

	qbvh->indices = malloc(size*sizeof(uint32_t));
	qbvh->indices_length = size;
	for (uint32_t i = 0; i != size; ++i)
		qbvh->indices[i] = i;

	qbvh->min_cells_per_leaf = min_cells_per_leaf;
	return qbvh;
}


// TODO do not copy bboxs and sort the original ones
void
bvh_build(struct Bvh *restrict qbvh, const float *restrict bboxs_orig)
{
	struct timespec tp0, tp1;
	// copy
	clock_gettime(CLOCK_MONOTONIC, &tp0);
	float *restrict bboxs;
	// 32 byte boundary alignment
	const size_t size = (4*qbvh->indices_length + 7) >> 3 << 3;
	posix_memalign(&bboxs, 32, sizeof(float[size]));
	memcpy(bboxs, bboxs_orig, sizeof(float[size]));
	clock_gettime(CLOCK_MONOTONIC, &tp1);
	printf("BVH: copy time %f ms\n", ((tp1.tv_sec*1e9 + tp1.tv_nsec) - (tp0.tv_sec*1e9 + tp0.tv_nsec))*1e-6);

	// TODO alignment on 32 byte boundary
	struct Task {
		float mid[BVH_WIDTH]; // mid points along which to split
		uint32_t offset[BVH_WIDTH];
		uint32_t length[BVH_WIDTH]; // 4 bits in child are not enough
		uint32_t id;
		enum Axis axis;
		uint8_t width; // current width of the node
	} stack[STACK_LENGTH];
	uint8_t stack_length = 0;

	// compute initial bbox
    clock_gettime(CLOCK_MONOTONIC, &tp0);
	float range[2] = {FLT_MAX, -FLT_MAX};
	for (uint32_t i = 0; i != qbvh->indices_length; ++i) {
		range[0] = MIN(range[0], bboxs[4*i + AXIS_X]);
		range[1] = MAX(range[1], bboxs[4*i + 2 + AXIS_X]);
	}
	clock_gettime(CLOCK_MONOTONIC, &tp1);
	printf("BVH: midpoint time %f ms\n", ((tp1.tv_sec*1e9 + tp1.tv_nsec) - (tp0.tv_sec*1e9 + tp0.tv_nsec))*1e-6);

	// initial task (root node)
	stack[0].mid[0] = 0.5f*(range[0] + range[1]);
    stack[0].offset[0] = 0;
    stack[0].length[0] = qbvh->indices_length;
	stack[0].id = qbvh->length++;
	stack[0].axis = AXIS_X;
	stack[0].width = 1;
    ++stack_length;

	/*
		example of array layout during construction of 8-way tree

		{n, __,  __, __} // 1 child (x-axis split)
		{n0, __, n1, __} // 2 children (y-axis split)
		{n00, n01, n10, n11} // 4 children (x-axis split)
		{n000, n001, n010, n011, n100, n101, n110, n111} // 8 children
	*/

	while (stack_length) {
		assert(stack_length - BVH_WIDTH < STACK_LENGTH);
		// TODO do not copy task around
		// BE EXTREMELY CAREFUL AS IF ITS NOT A COPY IT CAN GET OVERWRITTEN
		struct Task task = stack[--stack_length];

		// same for all children in the task
		const enum Axis axis = (task.axis + 1)%N_AXES;
		const uint8_t width = task.width;
		// for all children if less than with subdivide otherwise create
		for (int i = 0, step = BVH_WIDTH/width; i != BVH_WIDTH; i += step) {
			// split along midpoint
			float bbox[2][4];
			const uint32_t offset = split(bbox,
										  bboxs + 4*task.offset[i],
										  qbvh->indices + task.offset[i],
										  task.length[i],
										  task.axis,
										  task.mid[i]);
			const int padding = step/2;

			// compute left and right midpoints (arithmetic mean)
            task.mid[i] = 0.5f*(bbox[0][axis] + bbox[0][2 + axis]);
            task.mid[i + padding] = 0.5f*(bbox[1][axis] + bbox[1][2 + axis]);

            // FIRST COMPUTE THIS AS WE ARE OVERWRITTING THE LENGTH
            task.length[i + padding] = task.length[i] - offset;
            task.length[i] = offset;
            task.offset[i + padding] = task.offset[i] + offset;

            task.axis = axis;
            task.width = width*2;

			// padding == 1
			// reached width
			if (task.width == BVH_WIDTH) {
				// update left and right bbox
				qbvh->nodes[task.id].bbox[0][i] = bbox[0][0];
				qbvh->nodes[task.id].bbox[0][i + padding] = bbox[1][0];
				qbvh->nodes[task.id].bbox[1][i] = bbox[0][1];
				qbvh->nodes[task.id].bbox[1][i + padding] = bbox[1][1];
				qbvh->nodes[task.id].bbox[2][i] = bbox[0][2];
				qbvh->nodes[task.id].bbox[2][i + padding] = bbox[1][2];
				qbvh->nodes[task.id].bbox[3][i] = bbox[0][3];
				qbvh->nodes[task.id].bbox[3][i + padding] = bbox[1][3];

				if (qbvh->length + 2 > qbvh->capacity) {
            		qbvh->capacity *= 2;
		            qbvh->nodes = realloc(qbvh->nodes, qbvh->capacity*sizeof(struct BvhNode));
      			}

				// create left child
				if (task.length[i] > qbvh->min_cells_per_leaf) {
					qbvh->nodes[task.id].child[i] = qbvh->length;
					struct Task t;
					t.mid[0] = task.mid[i];
					t.offset[0] = task.offset[i];
					t.length[0] = task.length[i];
					t.axis = task.axis;
					t.width = 1;
					t.id = qbvh->length++;
					stack[stack_length++] = t;
				// finalize leaf
				} else {
					// subtract 1 from the length to use all 16 bits
					if (task.length[i])
						qbvh->nodes[task.id].child[i] = leaf_empty | task.length[i] - 1 << 27 | task.offset[i];
					else
						qbvh->nodes[task.id].child[i] = leaf_empty;
				}

				// create right child
				if (task.length[i + padding] > qbvh->min_cells_per_leaf) {
					qbvh->nodes[task.id].child[i + padding] = qbvh->length;
                    struct Task t;
                    t.mid[0] = task.mid[i + padding];
                    t.offset[0] = task.offset[i + padding];
                    t.length[0] = task.length[i + padding];
                    t.axis = task.axis;
                    t.width = 1;
                    t.id = qbvh->length++;
                    stack[stack_length++] = t;
                // finalize leaf
                } else {
					// subtract 1 from the length to use all 16 bits
                    if (task.length[i + padding])
                        qbvh->nodes[task.id].child[i + padding] = leaf_empty | task.length[i + padding] - 1 << 27 | task.offset[i + padding];
                    else
                        qbvh->nodes[task.id].child[i + padding] = leaf_empty;
				}
			}
		}

		// push back for processing
		if (task.width != BVH_WIDTH)
			stack[stack_length++] = task;
	}

	free(bboxs);
}


void
bvh_free(struct Bvh *qbvh)
{
	free(qbvh->nodes);
	free(qbvh->indices);
	free(qbvh);
}


uint32_t
bvh_intersect(const struct Bvh *restrict qbvh,
              uint32_t *restrict buffer,
              const uint32_t size,
              const float *restrict line)
{
	uint32_t total = 0;

	uint32_t stack[STACK_LENGTH];

	_Alignas(32) const float o[2] = {line[0], line[1]};
	_Alignas(32) const float d_inv[2] = {1.0f/(line[2] - o[0]), 1.0f/(line[3] - o[1])};

	// put root node onto the stack
	stack[0] = 0;
	uint8_t stack_length = 1;
	// check if type can hold the stack length
	assert(0x01 << 8*sizeof(stack_length) > STACK_LENGTH);

	_Alignas(32) int isect[BVH_WIDTH];

	while (stack_length) {
		assert(stack_length < STACK_LENGTH);
		const struct BvhNode *n = qbvh->nodes + stack[--stack_length];

#if BVH_WIDTH == 8
		intersect_8(isect, n->bbox, o, d_inv);
#else
		for (int i = 0; i != BVH_WIDTH; ++i)
			intersect(isect, n->bbox, o, d_inv, i);
#endif

		for (int i = 0; i != BVH_WIDTH; ++i) {
			if (n->child[i] == leaf_empty)
				continue;

			// line intersects or is inside
			if (isect[i]) {
				if (n->child[i] < 0) {
					assert(total < size);
					const uint8_t len = child_length(n->child[i]);

					// TODO use AVX for the copy
					for (int j = 0; j != len; ++j)
						buffer[total + j] = qbvh->indices[(n->child[i] & offset_mask) + j];
					total += len;
				} else
					stack[stack_length++] = n->child[i];
			}
		}
	}

	return total;
}


void
bvh_print(const struct Bvh *bvh)
{
	printf("BVH: nodes length %"PRIu32"\n\n", bvh->length);

	for (uint32_t i = 0; i != bvh->length; ++i) {
		printf("Node %"PRIu32"\n", i);
		const struct BvhNode *n = bvh->nodes + i;
		for (int j = 0; j != BVH_WIDTH; ++j) {
			printf("\tchild[%d] - leaf %d, offset %"PRIu32", "
			       "length %"PRIu8", bbox min(%f,%f) max(%f,%f)\n",
			       j,
			       (n->child[j] < 0),
			       n->child[j] & offset_mask,
			       child_length(n->child[j]),
			       n->bbox[0][j],
			       n->bbox[1][j],
			       n->bbox[2][j],
			       n->bbox[3][j]);
		}
		putchar('\n');
	}

}

void
bvh_stats(const struct Bvh *bvh)
{
	printf("BVH: width of the tree %d\n", BVH_WIDTH);
	printf("BVH: number of nodes %"PRIu32"\n", bvh->length);
	printf("BVH: number of indices %"PRIu32"\n", bvh->indices_length);
	printf("BVH: memory usage %zu bytes\n", (bvh->length*sizeof(struct BvhNode) +
	       bvh->indices_length*sizeof(uint32_t) + sizeof(struct Bvh)));
	printf("BVH: minimum leaf tet number %d\n", bvh->min_cells_per_leaf);

	uint32_t empty_count = 0;
	uint8_t max_depth = 0;
	struct Item {
		uint32_t id;
		uint8_t depth;
	} stack[STACK_LENGTH];
	stack[0].id = 0;
	stack[0].depth = 0;
	uint8_t stack_length = 1;

	while (stack_length) {
		assert(stack_length - BVH_WIDTH < STACK_LENGTH);
		const struct Item item = stack[--stack_length];
		for (int i = 0; i != BVH_WIDTH; ++i) {
			if (bvh->nodes[item.id].child[i] < 0) {
				max_depth = (item.depth > max_depth) ? item.depth : max_depth;
				empty_count += bvh->nodes[item.id].child[i] == leaf_empty;
			} else {
				stack[stack_length].id = bvh->nodes[item.id].child[i];
				stack[stack_length].depth = item.depth + 1;
				++stack_length;
			}
		}
	}

	printf("BVH: empty child count %"PRIu32"\n", empty_count);
	printf("BVH: maximum depth of the tree %"PRIu8"\n", max_depth);
}
