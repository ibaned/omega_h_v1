#ifndef UP_FROM_DOWN_H
#define UP_FROM_DOWN_H

struct up_from_down {
  unsigned* offsets;
  unsigned* edges;
  unsigned* directions;
};

struct up_from_down get_up_from_down(
    unsigned up_dim,
    unsigned down_dim,
    unsigned nup,
    unsigned ndown,
    unsigned* edges);

#endif
