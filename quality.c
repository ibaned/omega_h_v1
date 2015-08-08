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
