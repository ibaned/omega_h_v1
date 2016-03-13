#ifndef SUBSET_HPP
#define SUBSET_HPP

namespace omega_h {

struct mesh;

void tags_subset(struct mesh* in, struct mesh* out,
    unsigned dim, unsigned const* offsets);

struct mesh* subset_mesh(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* offsets);

void subset_verts_of_doms(
    struct mesh* m,
    unsigned dom_dim,
    unsigned const* offset_of_doms,
    unsigned* verts_of_prods);

}

#endif
