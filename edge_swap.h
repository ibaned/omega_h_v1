#ifndef EDGE_SWAP_H
#define EDGE_SWAP_H

#define MAX_EDGE_SWAP 7

extern unsigned const swap_mesh_sizes[MAX_EDGE_SWAP+1];
extern unsigned const swap_mesh_counts[MAX_EDGE_SWAP+1];

typedef unsigned const swap_tri_t[3];
extern swap_tri_t const* const swap_triangles[MAX_EDGE_SWAP+1];
extern unsigned const* const swap_meshes[MAX_EDGE_SWAP+1];
extern unsigned const* const* const swap_int_edges[MAX_EDGE_SWAP + 1];
extern unsigned const swap_nint_edges[MAX_EDGE_SWAP + 1];

struct swap_choice {
  unsigned code;
  int padding__;
  double quality;
};

struct swap_choice choose_edge_swap(
    unsigned ring_size,
    double (*edge_x)[3],
    double (*ring_x)[3]);

void get_swap_tets(
    unsigned ring_size,
    unsigned code,
    unsigned const edge_v[2],
    unsigned const* ring_v,
    unsigned* out);

void get_swap_edges(
    unsigned ring_size,
    unsigned code,
    unsigned const* ring_v,
    unsigned* out);

void get_swap_tris(
    unsigned ring_size,
    unsigned code,
    unsigned const edge_v[2],
    unsigned const* ring_v,
    unsigned* out);

#endif
