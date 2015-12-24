#ifndef CONCAT_H
#define CONCAT_H

void concat_verts_of_ents(
    unsigned ent_dim,
    unsigned nents,
    unsigned ngen_ents,
    unsigned const* verts_of_ents,
    unsigned const* offset_of_same_ents,
    unsigned const* verts_of_gen_ents,
    unsigned* nents_out,
    unsigned** verts_of_ents_out);

#endif
