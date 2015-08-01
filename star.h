#ifndef STAR_H
#define STAR_H

#include "up_from_down.h"

struct star {
  unsigned* offsets;
  unsigned* edges;
};

struct star get_star(
    unsigned down_dim,
    unsigned up_dim,
    unsigned ndown,
    unsigned const* up_offsets,
    unsigned const* up_edges,
    unsigned const* down_edges);

#endif
