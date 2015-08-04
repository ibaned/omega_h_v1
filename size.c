#include "size.h"
#include "algebra.h"

double edge_length(double const coords[2][3])
{
  return vector_distance(coords[1], coords[0], 3);
}

double triangle_area(double const coords[3][3])
{
  double v[2][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return vector_norm(x, 3) / 2.0;
}

double tet_volume(double const coords[4][3])
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
