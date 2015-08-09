#ifndef COARSEN_REDUCED_H
#define COARSEN_REDUCED_H

int coarsen_reduced(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    unsigned** p_class_dim,
    double (*size_function)(double const x[]));

#endif
