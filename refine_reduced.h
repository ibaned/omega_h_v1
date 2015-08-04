#ifndef REFINE_REDUCED_H
#define REFINE_REDUCED_H

struct reduced_mesh {
  unsigned nelem;
  unsigned nvert;
  unsigned* elem_verts;
  double* coords;
};

struct reduced_mesh refine_reduced(
    unsigned elem_dim,
    unsigned nelem,
    unsigned nvert,
    unsigned const* elem_verts,
    double const* coords,
    double (*size_function)(double const x[]));

#endif
