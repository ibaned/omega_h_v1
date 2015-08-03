#include "measure_edges.h"
#include "algebra.h"
#include "size.h"
#include <stdlib.h>

double* measure_edges(
    unsigned nedges,
    unsigned const* edge_verts,
    double const* coords,
    double const* size)
{
  double* out = malloc(sizeof(double) * nedges);
  for (unsigned i = 0; i < nedges; ++i) {
    unsigned const* edge_vert = edge_verts + i * 2;
    double edge_coord[2][3];
    copy_vector(coords + edge_vert[0] * 3, edge_coord[0], 3);
    copy_vector(coords + edge_vert[1] * 3, edge_coord[1], 3);
    double length = edge_length(edge_coord);
    double desired_length = (size[edge_vert[0]] + size[edge_vert[1]]) / 2;
    out[i] = length / desired_length;
  }
  return out;
}
