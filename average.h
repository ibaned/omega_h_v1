#ifndef AVERAGE_H
#define AVERAGE_H

static inline void average_element_field(
    unsigned nelem_vert,
    unsigned const* elem_vert,
    double const* field,
    unsigned ncomp,
    double* out)
{
  for (unsigned i = 0; i < ncomp; ++i)
    out[i] = 0;
  for (unsigned i = 0; i < nelem_vert; ++i) {
    unsigned vert = elem_vert[i];
    for (unsigned j = 0; j < ncomp; ++j)
      out[j] += field[vert * ncomp + j];
  }
  for (unsigned i = 0; i < ncomp; ++i)
    out[i] /= nelem_vert;
}

#endif
