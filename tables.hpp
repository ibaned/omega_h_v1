#ifndef TABLES_H
#define TABLES_H

#include "loop.hpp"

#define MAX_DOWN 6
#define MAX_UP 128
#define INVALID (~((unsigned)0))

extern unsigned const* const the_box_conns[4];
extern double const* const the_box_coords[4];
extern unsigned const the_box_nelems[4];
extern unsigned const the_box_nverts[4];
extern unsigned const the_down_degrees[4][4];
LOOP_CONST extern unsigned const* const* const* const* const the_canonical_orders[4];
LOOP_CONST extern unsigned const* const* const the_opposite_orders[4];

static inline unsigned get_opposite_dim(
    unsigned elem_dim,
    unsigned ent_dim)
{
  return elem_dim - 1 - ent_dim;
}

char const* get_ent_name(unsigned dim, unsigned long n);

#endif
