#ifndef COARSEN_BY_SIZE_H
#define COARSEN_BY_SIZE_H

unsigned coarsen_by_size(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    unsigned** p_class_dim,
    double (*size_function)(double const x[]),
    double quality_floor);

#endif
