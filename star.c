#include "star.h"
#include "tables.h"
#include "ints.h"
#include <stdlib.h>

static unsigned add_unique(unsigned* a, unsigned n, unsigned e)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] == e)
      return n;
  a[n] = e;
  return n + 1;
}

static unsigned get_ent_star(
    unsigned const* up_offsets,
    unsigned const* up_edges,
    unsigned const* down_edges,
    unsigned down_degree,
    unsigned down,
    unsigned* ent_star)
{
  unsigned first_up = up_offsets[down];
  unsigned last_up = up_offsets[down + 1];
  unsigned size = 0;
  for (unsigned i = first_up; i < last_up; ++i) {
    unsigned up = up_edges[i];
    for (unsigned j = 0; j < down_degree; ++j) {
      unsigned star_down = down_edges[up * down_degree + j];
      size = add_unique(ent_star, size, star_down);
    }
  }
  return size;
}

struct star get_star(
    unsigned down_dim,
    unsigned up_dim,
    unsigned ndown,
    unsigned const* up_offsets,
    unsigned const* up_edges,
    unsigned const* down_edges)
{
  struct star out;
  unsigned star_buf[MAX_UP * MAX_DOWN];
  unsigned down_degree = the_down_degrees[up_dim][down_dim];
  unsigned* degrees = malloc(sizeof(unsigned) * ndown);
  ints_zero(degrees, ndown);
  for (unsigned i = 0; i < ndown; ++i)
    degrees[i] = get_ent_star(
        up_offsets,
        up_edges,
        down_edges,
        down_degree,
        i,
        star_buf);
  out.offsets = ints_exscan(degrees, ndown);
  unsigned nedges = out.offsets[ndown];
  out.edges = malloc(sizeof(unsigned) * nedges);
  for (unsigned i = 0; i < ndown; ++i) {
    get_ent_star(
        up_offsets,
        up_edges,
        down_edges,
        down_degree,
        i,
        star_buf);
    unsigned first_star = out.offsets[i];
    unsigned last_star = out.offsets[i + 1];
    for (unsigned j = first_star; j < last_star; ++j)
      out.edges[j] = star_buf[j - first_star];
  }
  return out;
}
