#ifndef COLLAPSE_CLASSIF_H
#define COLLAPSE_CLASSIF_H

void check_collapse_classif(
    unsigned elem_dim,
    unsigned nedges,
    unsigned* col_codes,
    unsigned const* class_dim_of_verts,
    unsigned const* verts_of_elems,
    unsigned const* verts_of_edges,
    unsigned const* verts_of_verts_offsets,
    unsigned const* verts_of_verts,
    unsigned const* elems_of_edges_offsets,
    unsigned const* elems_of_edges,
    unsigned const* elems_of_edges_directions);

#endif
