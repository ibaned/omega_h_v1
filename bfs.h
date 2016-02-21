#ifndef BFS_H
#define BFS_H

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

#endif
