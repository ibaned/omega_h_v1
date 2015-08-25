#include "inside.h"
#include "algebra.h"
#include "jacobian.h"

void point_to_edge(double edge[2][3], double const pt[3], double b[2])
{
  double v[2][3];
  subtract_vectors(edge[1], edge[0], v[0], 3);
  subtract_vectors(pt, edge[0], v[1], 3);
  b[0] = dot_product(v[1], v[0], 3) / dot_product(v[0], v[0], 3);
  b[1] = 1.0 - b[0];
}

void point_to_triangle(double tri[3][3], double const pt[3], double b[3])
{
  double jact[3][3];
  subtract_vectors(tri[1], tri[0], jact[0], 3);
  subtract_vectors(tri[2], tri[0], jact[1], 3);
  cross_product(jact[0], jact[1], jact[2]);
  double jac[3][3];
  transp_3x3(jact, jac);
  double jaci[3][3];
  invert_3x3(jac, jaci);
  double v[3];
  subtract_vectors(pt, tri[0], v, 3);
  mv_3x3(jaci, v, b);
  b[2] = 1.0 - b[0] - b[1];
}

void point_to_tet(double tet[4][3], double const pt[3], double b[4])
{
  double jac[3][3];
  element_jacobian(3, tet, jac);
  double jaci[3][3];
  invert_3x3(jac, jaci);
  double v[3];
  subtract_vectors(pt, tet[0], v, 3);
  mv_3x3(jaci, v, b);
  b[3] = 1.0 - b[0] - b[1] - b[2];
}

inside_function const the_inside_functions[4] = {
  0,
  point_to_edge,
  point_to_triangle,
  point_to_tet
};
