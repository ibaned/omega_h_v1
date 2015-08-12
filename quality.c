#include "quality.h"
#include "size.h"
#include "algebra.h"
#include "tables.h"

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
  double a = triangle_area(coords);
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


enum tri_qual triangle_quality_type(
    double coords[3][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_edge_out)
{
  double q = triangle_quality(coords);
  if (q >= qual_floor)
    return GOOD_TRI;
  unsigned const* const* fev = the_canonical_orders[2][1][0];
  double l[3];
  for (unsigned i = 0; i < 3; ++i)
    l[i] = vector_distance(
        coords[fev[i][1]], coords[fev[i][0]], 3);
  double maxl = l[0];
  double minl = l[0];
  for (unsigned i = 1; i < 3; ++i) {
    if (l[i] > maxl)
      maxl = l[i];
    if (l[i] < minl)
      minl = l[i];
  }
  double r = minl / maxl;
  if (r < edge_ratio_floor)
    return SHORT_EDGE_TRI;
  if (!key_edge_out)
    return SLIVER_TRI;
  double v[2][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  /* scalar projection of (v2 - v0) onto (v1 - v0) */
  double sp = dot_product(v[1], v[0], 3) / l[0];
  if (sp < 0)
    *key_edge_out = 1;
  else if (sp > l[0])
    *key_edge_out = 2;
  else
    *key_edge_out = 0;
  return SLIVER_TRI;
}

