#include "quality.h"
#include "size.h"
#include "algebra.h"
#include "tables.h"
#include "inside.h"
#include "doubles.h"
#include "mesh.h"

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
  double a = triangle_area(coords);
  double quality = a / (lrms * lrms);
  return quality / PERFECT_TRIANGLE_QUALITY;
}

double triangle_xy_quality(double coords[3][3])
{
  unsigned const* const* fev = the_canonical_orders[2][1][0];
  double sum_lsq = 0;
  for (unsigned i = 0; i < 3; ++i) {
    double lsq = vector_squared_distance(
        coords[fev[i][1]], coords[fev[i][0]], 3);
    sum_lsq += lsq;
  }
  double lrms = sqrt(sum_lsq / 3);
  /* use the area of the triangle projected onto the XY plane,
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
  double arms = sqrt(sum_asq / 4.);
  double v = tet_volume(coords);
  double root_arms = sqrt(arms);
  double quality = v / CUBE(root_arms);
  return quality / PERFECT_TET_QUALITY;
}

static double one(double x[][3])
{
  (void) x;
  return 1;
}

quality_function const the_quality_functions[4] = {
  one,
  one,
  triangle_quality,
  tet_quality
};

quality_function const the_equal_order_quality_functions[4] = {
  one,
  one,
  triangle_xy_quality,
  tet_quality
};

double* element_qualities(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords)
{
  double* out = malloc(sizeof(double) * nelems);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  quality_function qf = the_equal_order_quality_functions[elem_dim];
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

double* mesh_qualities(struct mesh* m)
{
  return element_qualities(mesh_dim(m),
      mesh_count(m, mesh_dim(m)),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_nodal_field(m, "coordinates")->data);
}

double mesh_min_quality(struct mesh* m)
{
  return min_element_quality(mesh_dim(m),
      mesh_count(m, mesh_dim(m)),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_nodal_field(m, "coordinates")->data);
}
