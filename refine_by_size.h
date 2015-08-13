#ifndef REFINE_BY_SIZE_H
#define REFINE_BY_SIZE_H

unsigned refine_by_size(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    double (*size_function)(double const x[]));

#endif
