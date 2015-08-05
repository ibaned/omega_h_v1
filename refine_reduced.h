#ifndef REFINE_REDUCED_H
#define REFINE_REDUCED_H

void refine_reduced(
    unsigned elem_dim,
    unsigned nelem,
    unsigned nvert,
    unsigned const* elem_verts,
    double const* coords,
    double (*size_function)(double const x[]),
    unsigned* nelem_out,
    unsigned* nverts_out,
    unsigned** verts_of_elems_out,
    double** coords_out);

#endif
