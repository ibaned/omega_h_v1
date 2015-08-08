#ifndef COLLAPSES_TO_VERTS_H
#define COLLAPSES_TO_VERTS_H

void collapses_to_verts(
    unsigned nverts,
    unsigned const* verts_of_edges,
    unsigned const* edges_of_verts_offsets,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_directions,
    unsigned const* col_codes,
    double const* col_quals_of_edges,
    unsigned** gen_offset_of_verts_out,
    unsigned** gen_vert_of_verts_out,
    double** col_qual_of_verts_out);

#endif
