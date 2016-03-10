#ifndef SIZE_H
#define SIZE_H

#include "algebra.h"

LOOP_INOUT static inline double
edge_length(double coords[2][3])
{
  return vector_distance(coords[1], coords[0], 3);
}

LOOP_INOUT static inline double
triangle_area(double coords[3][3])
{
  double v[2][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return vector_norm(x, 3) / 2.0;
}

LOOP_INOUT static inline double
triangle_z_area(double coords[3][3])
{
  double v[2][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return x[2] / 2.0;
}

LOOP_INOUT static inline double
tet_volume(double coords[4][3])
{
  double v[3][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  subtract_vectors(coords[3], coords[0], v[2], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return dot_product(x, v[2], 3) / 6.0;
}

LOOP_INOUT static inline double
measure_entity(unsigned dim, double (*coords)[3])
{
  switch (dim) {
    case 3: return tet_volume(coords);
    case 2: return triangle_area(coords);
    case 1: return edge_length(coords);
  }
  return 0.0;
}

double* element_sizes(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords);

struct mesh;

double* mesh_element_sizes(struct mesh* m);

double* mesh_measure_edges_for_adapt(struct mesh* m);

void mesh_identity_size_field(struct mesh* m, char const* output_name);

#endif
