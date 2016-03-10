#ifndef BFS_HPP
#define BFS_HPP

/* all arrays on HOST ! */

void bfs_continue(
    unsigned* queue,
    unsigned* p_begin,
    unsigned* p_end,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned* comp,
    unsigned* layer);

void bfs_from(
    unsigned start,
    unsigned the_comp,
    unsigned* queue,
    unsigned* p_begin,
    unsigned* p_end,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned* comp,
    unsigned* layer);

void bfs_full(
    unsigned n,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned* comp,
    unsigned* layer,
    unsigned* sorted);

void connected_components(
    unsigned n,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned* comp);

#endif
