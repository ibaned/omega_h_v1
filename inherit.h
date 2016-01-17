#ifndef INHERIT_H
#define INHERIT_H

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

void concat_verts_of_elems(
    unsigned elem_dim,
    unsigned nelems,
    unsigned ngen_elems,
    unsigned const* verts_of_elems,
    unsigned const* offset_of_same_elems,
    unsigned* verts_of_gen_elems,
    unsigned* nelems_out,
    unsigned** verts_of_elems_out);

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
