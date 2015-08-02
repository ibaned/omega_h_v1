#ifndef REFINE_TOPOLOGY_H
#define REFINE_TOPOLOGY_H

struct refined_topology {
  unsigned nent;
  unsigned padding__;
  unsigned* ent_verts;
};

struct refined_topology refine_topology(
    unsigned elem_dim,
    unsigned split_dim,
    unsigned ent_dim,
    unsigned nelem,
    unsigned const* elem_verts,
    unsigned const* elem_split_offset,
    unsigned const* elem_split_vert,
    unsigned const* elem_split_direction);

#endif
