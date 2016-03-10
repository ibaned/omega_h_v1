#ifndef AVERAGE_H
#define AVERAGE_H

#include "loop.hpp"

LOOP_INOUT static inline void
average_element_field(
    unsigned verts_per_elem,
    unsigned const* verts_of_elem,
    double const* field,
    unsigned comps_per_vert,
    double* out)
{
  for (unsigned i = 0; i < comps_per_vert; ++i)
    out[i] = 0;
  for (unsigned i = 0; i < verts_per_elem; ++i) {
    unsigned vert = verts_of_elem[i];
    for (unsigned j = 0; j < comps_per_vert; ++j)
      out[j] += field[vert * comps_per_vert + j];
  }
  for (unsigned i = 0; i < comps_per_vert; ++i)
    out[i] /= verts_per_elem;
}

double* interp_to_elems(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned ncomps,
    double const* in);

struct mesh;

void mesh_interp_to_elems(struct mesh* m, char const* name);

#endif
