#ifndef COARSEN_QUALITIES
#define COARSEN_QUALITIES

double* coarsen_qualities(
    unsigned elem_dim,
    unsigned nedges,
    unsigned* col_codes,
    unsigned const* verts_of_elems,
    unsigned const* verts_of_edges,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    unsigned const* elems_of_verts_directions,
    double const* coords,
    double quality_limit);

#endif
