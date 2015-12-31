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

void inherit_uint_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned prod_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4],
    char const* name);

void refine_class(
    struct mesh* m_in,
    struct mesh* m_out,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

#endif
