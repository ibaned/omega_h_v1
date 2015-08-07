#ifndef REFINE_REDUCED_H
#define REFINE_REDUCED_H

int refine_reduced(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    double (*size_function)(double const x[]));

#endif
