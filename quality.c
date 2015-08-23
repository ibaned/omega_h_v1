#include "quality.h"
#include "size.h"
#include "algebra.h"
#include "tables.h"
#include "inside.h"
#include "doubles.h"

#include <assert.h>
#include <stdlib.h>

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

/* These normalization constants are the values of
 * the original quality measures for regular simplices.
 * The formulae for area and volume of regular simplices
 * based on edge length are on the Wikipedia.
 */
#define PERFECT_TRIANGLE_QUALITY (sqrt(3.0) / 4.0)
#define CUBE(x) ((x)*(x)*(x))
#define PERFECT_TET_QUALITY \
  ((sqrt(2.0) / 12.0) / CUBE(sqrt(PERFECT_TRIANGLE_QUALITY)))

double triangle_quality(double coords[3][3])
{
  unsigned const* const* fev = the_canonical_orders[2][1][0];
  double sum_lsq = 0;
  for (unsigned i = 0; i < 3; ++i) {
    double lsq = vector_squared_distance(
        coords[fev[i][1]], coords[fev[i][0]], 3);
    sum_lsq += lsq;
  }
  double lrms = sqrt(sum_lsq / 3);
  /* highly dubious.
     I yelled at someone for doing this and now I'm doing it.
     use the area of the triangle projected onto the XY plane,
     with the sign of the normal's Z component.
     the point of this is the prevent flipping triangles */
  double a = triangle_z_area(coords);
  double quality = a / (lrms * lrms);
  return quality / PERFECT_TRIANGLE_QUALITY;
}

double tet_quality(double coords[4][3])
{
  unsigned const* const* rfv = the_canonical_orders[3][2][0];
  double sum_asq = 0;
  for (unsigned i = 0; i < 4; ++i) {
    double tri_coords[3][3];
    for (unsigned j = 0; j < 3; ++j)
      copy_vector(coords[rfv[i][j]], tri_coords[j], 3);
    double a = triangle_area(tri_coords);
    sum_asq += a * a;
  }
  double arms = sqrt(sum_asq / 4);
  double v = tet_volume(coords);
  double root_arms = sqrt(arms);
  double quality = v / CUBE(root_arms);
  return quality / PERFECT_TET_QUALITY;
}

quality_function const the_quality_functions[4] = {
  0,
  0,
  triangle_quality,
  tet_quality
};

static void edge_lengths(
    unsigned elem_dim,
    double coords[][3],
    double l[])
{
  unsigned ne = the_down_degrees[elem_dim][1];
  unsigned const* const* lev = the_canonical_orders[elem_dim][1][0];
  for (unsigned i = 0; i < ne; ++i)
    l[i] = vector_distance(
        coords[lev[i][1]], coords[lev[i][0]], 3);
}

static double edge_length_ratio(
    double l[],
    unsigned n)
{
  double maxl = l[0];
  double minl = l[0];
  for (unsigned i = 1; i < 3; ++i) {
    if (l[i] > maxl)
      maxl = l[i];
    if (l[i] < minl)
      minl = l[i];
  }
  return minl / maxl;
}

enum quality_type triangle_quality_type(
    double coords[3][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_edge_out)
{
  double q = triangle_quality(coords);
  if (q >= qual_floor)
    return GOOD_ELEM;
  double l[3];
  edge_lengths(2, coords, l);
  double r = edge_length_ratio(l, 3);
  if (r < edge_ratio_floor)
    return SHORT_EDGE_ELEM;
  if (!key_edge_out)
    return SLIVER_ELEM;
  /* project the cap vertex onto the base edge */
  /* the first two vertices of the triangle are
     ordered the same as the vertices of its first edge */
  double b[2];
  point_to_edge(coords, coords[2], b);
  if (b[0] < 0)
    *key_edge_out = 1;
  else if (b[0] > 1)
    *key_edge_out = 2;
  else
    *key_edge_out = 0;
  return SLIVER_ELEM;
}

static struct {
  enum quality_type qt;
  unsigned ke;
} const tet_quality_table[8] = {
  {0,0}, /* impossible for positive volume */
  {CAP_TET,     1}, /* 001 */
  {CAP_TET,     2}, /* 010 */
  {SLIVER_ELEM, 1}, /* 011 */
  {CAP_TET,     3}, /* 100 */
  {SLIVER_ELEM, 0}, /* 101 */
  {SLIVER_ELEM, 2}, /* 110 */
  {CAP_TET,     0}, /* 111 */
};

enum quality_type tet_quality_type(
    double coords[4][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_ent_out)
{
  double q = tet_quality(coords);
  if (q >= qual_floor)
    return GOOD_ELEM;
  double l[6];
  edge_lengths(3, coords, l);
  double r = edge_length_ratio(l, 6);
  if (r < edge_ratio_floor)
    return SHORT_EDGE_ELEM;
  /* project the cap vertex onto the base triangle */
  /* the first three vertices of the tet are
     ordered the same as the vertices of its first triangle */
  double b[3];
  point_to_triangle(coords, coords[3], b);
  int c = 0;
  for (unsigned i = 0; i < 3; ++i)
    if (b[i] > 0)
      c |= (1<<i);
  assert(c);
  if (key_ent_out)
    *key_ent_out = tet_quality_table[c].ke;
  return tet_quality_table[c].qt;
}

quality_type_function const the_quality_type_functions[4] = {
  0,
  0,
  triangle_quality_type,
  tet_quality_type
};

double* element_qualities(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords)
{
  double* out = malloc(sizeof(double) * nelems);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  quality_function qf = the_quality_functions[elem_dim];
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    double elem_x[MAX_DOWN][3];
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned vert = verts_of_elem[j];
      copy_vector(coords + vert * 3, elem_x[j], 3);
    }
    out[i] = qf(elem_x);
  }
  return out;
}

double min_element_quality(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords)
{
  double* quals = element_qualities(elem_dim, nelems, verts_of_elems, coords);
  double mq = doubles_min(quals, nelems);
  free(quals);
  return mq;
}
