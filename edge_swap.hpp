#ifndef EDGE_SWAP_H
#define EDGE_SWAP_H

#include "loop.h"

#define MAX_EDGE_SWAP 7

LOOP_CONST extern unsigned const swap_mesh_sizes[MAX_EDGE_SWAP+1];
LOOP_CONST extern unsigned const swap_mesh_counts[MAX_EDGE_SWAP+1];

typedef unsigned const swap_tri_t[3];
LOOP_CONST extern swap_tri_t const* const swap_triangles[MAX_EDGE_SWAP+1];
LOOP_CONST extern unsigned const* const swap_meshes[MAX_EDGE_SWAP+1];

LOOP_CONST extern unsigned const* const* const swap_int_edges[MAX_EDGE_SWAP + 1];
LOOP_CONST extern unsigned const swap_nint_edges[MAX_EDGE_SWAP + 1];

struct swap_choice {
  unsigned code;
  int padding__;
  double quality;
};

LOOP_IN struct swap_choice choose_edge_swap(
    unsigned ring_size,
    double (*edge_x)[3],
    double (*ring_x)[3]);

LOOP_IN void get_swap_ents(
    unsigned ring_size,
    unsigned code,
    unsigned ent_dim,
    unsigned const edge_v[2],
    unsigned const* ring_v,
    unsigned* out);

LOOP_IN unsigned count_swap_ents(
    unsigned ring_size,
    unsigned ent_dim);

#endif
