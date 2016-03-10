#ifndef INHERIT_HPP
#define INHERIT_HPP

struct mesh;

void setup_refine(
    struct mesh* m,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned* gen_dom_offsets[4],
    /* out: */
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4],
    unsigned ngen_offsets[5]);

void setup_coarsen(
    struct mesh* m,
    unsigned ent_dim,
    unsigned* gen_offset_of_ents,
    unsigned* offset_of_same_ents,
    /* out: */
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

void setup_swap(
    struct mesh* m,
    unsigned ent_dim,
    unsigned* gen_offset_of_edges,
    unsigned* offset_of_same_ents,
    /* out: */
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

void inherit_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned prod_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

void concat_verts_of_ents(
    unsigned ent_dim,
    unsigned nents,
    unsigned ngen_ents,
    unsigned const* verts_of_ents,
    unsigned const* offset_of_same_ents,
    unsigned const* verts_of_gen_ents,
    unsigned* p_nents_out,
    unsigned** p_verts_of_ents_out);

unsigned* concat_uints_inherited(
    unsigned width,
    unsigned const ngen_offsets[5],
    unsigned* gen_data[4]);

double* concat_doubles_inherited(
    unsigned width,
    unsigned const ngen_offsets[5],
    double* gen_data[4]);

void make_ngen_from_doms(
    unsigned const ndoms[4],
    unsigned* prods_of_doms_offsets[4],
    unsigned ngen[4]);

void make_ngen_offsets(
    unsigned const ngen[4],
    /* out: */
    unsigned ngen_offsets[5]);

#endif
