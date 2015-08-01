#ifndef UP_FROM_DOWN_H
#define UP_FROM_DOWN_H

struct up_adj {
  unsigned* offsets;
  unsigned* edges;
  unsigned* directions;
};

struct up_adj up_from_down(
    unsigned up_dim,
    unsigned down_dim,
    unsigned nup,
    unsigned ndown,
    unsigned const* down_edges);

#endif
