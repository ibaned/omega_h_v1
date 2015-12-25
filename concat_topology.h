#ifndef CONCAT_H
#define CONCAT_H

void concat_verts_of_ents(
    unsigned ent_dim,
    unsigned nents,
    unsigned const* verts_of_ents,
    unsigned const* offset_of_same_ents,
    unsigned const ngen_ents[4],
    unsigned* verts_of_gen_ents[4],
    unsigned* p_nents_out,
    unsigned** p_verts_of_ents_out);

void concat_verts_of_elems(
    unsigned elem_dim,
    unsigned nelems,
    unsigned ngen_elems,
    unsigned const* verts_of_elems,
    unsigned const* offset_of_same_elems,
    unsigned* verts_of_gen_elems,
    unsigned* nelems_out,
    unsigned** verts_of_elems_out);

#endif
