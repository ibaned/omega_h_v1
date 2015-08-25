#include "size.h"
#include "algebra.h"
#include <stdlib.h>
#include <assert.h>

double edge_length(double coords[2][3])
{
  return vector_distance(coords[1], coords[0], 3);
}

double triangle_area(double coords[3][3])
{
  double v[2][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return vector_norm(x, 3) / 2.0;
}

double triangle_z_area(double coords[3][3])
{
  double v[2][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return x[2] / 2.0;
}

double tet_volume(double coords[4][3])
{
  double v[3][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  subtract_vectors(coords[3], coords[0], v[2], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return dot_product(x, v[2], 3) / 6.0;
}

element_measure const the_element_measures[4] = {
  0,
  edge_length,
  triangle_area,
  tet_volume
};

double* identity_size_field(
    unsigned nverts,
    unsigned const* vert_of_verts_offsets,
    unsigned const* vert_of_verts,
    double const* coords)
{
  double* out = malloc(sizeof(double) * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_use = vert_of_verts_offsets[i];
    unsigned end_use = vert_of_verts_offsets[i + 1];
    double edge_x[2][3];
    copy_vector(coords + i * 3, edge_x[0], 3);
    double max = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned ov = vert_of_verts[j];
      copy_vector(coords + ov * 3, edge_x[1], 3);
      double l = edge_length(edge_x);
      if (l > max)
        max = l;
    }
    out[i] = max;
  }
  return out;
}
