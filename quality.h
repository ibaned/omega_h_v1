#ifndef QUALITY_H
#define QUALITY_H

#include <assert.h>

#include "size.h"

/* These normalization constants are the values of
 * the original quality measures for regular simplices.
 * The formulae for area and volume of regular simplices
 * based on edge length are on the Wikipedia.
 */
#define PERFECT_TRIANGLE_QUALITY (sqrt(3.0) / 4.0)
#define CUBE(x) ((x)*(x)*(x))
#define PERFECT_TET_QUALITY \
  ((sqrt(2.0) / 12.0) / CUBE(sqrt(PERFECT_TRIANGLE_QUALITY)))

/* scale-invariant stiffness matrix conditioning quality measures
 * from:
 *
 * Shewchuk, J.
 * "What is a good linear finite element?
 *  interpolation, conditioning, anisotropy,
 *  and quality measures (preprint)."
 * University of California at Berkeley 73 (2002).
 *
 * Even when edges are normal length, these measures
 * are zero for zero volume and negative for inverted
 * elements.
 */

LOOP_INOUT static inline double
triangle_quality(double coords[3][3])
{
  double sum_lsq = 0;
  for (unsigned i = 0; i < 3; ++i) {
    double lsq = vector_squared_distance(
        coords[i], coords[(i + 1) % 3], 3);
    sum_lsq += lsq;
  }
  double lrms = sqrt(sum_lsq / 3);
  double a = triangle_area(coords);
  double quality = a / (lrms * lrms);
  return quality / PERFECT_TRIANGLE_QUALITY;
}

/* same as above, but for the triangle projected to
   the X-Y plane. this measure is negative if the
   normal points in the negative Z direction, which
   is used to prevent "negative" triangles in 2D space. */
LOOP_INOUT static inline double
triangle_xy_quality(double coords[3][3])
{
  double sum_lsq = 0;
  for (unsigned i = 0; i < 3; ++i) {
    double lsq = vector_squared_distance(
        coords[i], coords[(i + 1) % 3], 3);
    sum_lsq += lsq;
  }
  double lrms = sqrt(sum_lsq / 3);
  double a = triangle_z_area(coords);
  double quality = a / (lrms * lrms);
  return quality / PERFECT_TRIANGLE_QUALITY;
}

LOOP_INOUT static inline double
tet_quality(double coords[4][3])
{
  unsigned const rfv[4][3] = {
    {0,1,2},
    {1,0,3},
    {2,1,3},
    {0,2,3}};
  double sum_asq = 0;
  for (unsigned i = 0; i < 4; ++i) {
    double tri_coords[3][3];
    for (unsigned j = 0; j < 3; ++j)
      copy_vector(coords[rfv[i][j]], tri_coords[j], 3);
    double a = triangle_area(tri_coords);
    sum_asq += a * a;
  }
  double arms = sqrt(sum_asq / 4.);
  double v = tet_volume(coords);
  double root_arms = sqrt(arms);
  double quality = v / CUBE(root_arms);
  return quality / PERFECT_TET_QUALITY;
}

LOOP_INOUT static inline double
entity_quality(unsigned dim, double (*coords)[3])
{
  switch (dim) {
    case 3: return tet_quality(coords);
    case 2: return triangle_quality(coords);
  }
  assert(0);
#ifdef __CUDACC__
  return 0;
#endif
}

LOOP_INOUT static inline double
element_quality(unsigned dim, double (*coords)[3])
{
  switch (dim) {
    case 3: return tet_quality(coords);
    case 2: return triangle_xy_quality(coords);
  }
  assert(0);
#ifdef __CUDACC__
  return 0;
#endif
}

double* element_qualities(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords);

double min_element_quality(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords);

struct mesh;

double* mesh_qualities(struct mesh* m);
double mesh_min_quality(struct mesh* m);

#endif
