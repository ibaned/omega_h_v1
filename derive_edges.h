#ifndef DERIVE_EDGES_H
#define DERIVE_EDGES_H

struct derived_edges {
  unsigned nedge;
  unsigned padding__;
  unsigned* edge_verts;
};

struct derived_edges derive_edges(
    unsigned elem_dim,
    unsigned nelem,
    unsigned nvert,
    unsigned const* elem_verts);

#endif
