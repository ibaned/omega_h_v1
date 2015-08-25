#ifndef AVERAGE_H
#define AVERAGE_H

static inline void average_element_field(
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

#endif
