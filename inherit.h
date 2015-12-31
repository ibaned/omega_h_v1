#ifndef MODIFY_H
#define MODIFY_H

struct mesh;

void setup_refine(
    struct mesh* m,
    unsigned prod_dim,
    unsigned src_dim,
    unsigned* gen_dom_offsets[4],
    /* out: */
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4]);

void inherit_uint_tag(
    struct mesh* m,
    struct mesh* m_out,
    unsigned prod_dim,
    unsigned ndoms[4],
    unsigned* prods_of_doms_offsets[4],
    char const* name);

#endif
