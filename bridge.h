#ifndef BRIDGE_H
#define BRIDGE_H

struct bridged_graph {
  unsigned nedge;
  unsigned padding__;
  unsigned* edge_verts;
};

struct bridged_graph bridge_graph(
    unsigned nverts,
    unsigned const* offsets,
    unsigned const* edges);

#endif
