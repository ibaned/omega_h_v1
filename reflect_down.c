#include "reflect_down.h"

#include <assert.h>
#include <stdio.h>
#include <mpi.h>

#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

LOOP_INOUT static unsigned copy(
    unsigned const* a,
    unsigned* b,
    unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return n;
}

LOOP_INOUT static unsigned copy_except(
    unsigned const* a,
    unsigned* b,
    unsigned n,
    unsigned exclude_this)
{
  unsigned j = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] != exclude_this)
      b[j++] = a[i];
  return j;
}

LOOP_INOUT static unsigned intersect(
    unsigned* a,
    unsigned na,
    unsigned const* b,
    unsigned nb)
{
  unsigned j = 0;
  for (unsigned i = 0; i < na; ++i)
    if (has(b, nb, a[i]))
      a[j++] = a[i];
  return j;
}

LOOP_KERNEL(reflect_down_entity,
    unsigned dual_mode,
    unsigned const* verts_of_highs,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned verts_per_high,
    unsigned lows_per_high,
    unsigned verts_per_low,
    unsigned* lows_of_highs,
    unsigned const* const* high_verts_of_lows)

  unsigned const* verts_of_high = verts_of_highs + i * verts_per_high;
  unsigned* lows_of_high = lows_of_highs + i * lows_per_high;
  for (unsigned j = 0; j < lows_per_high; ++j) {
    unsigned const* high_verts_of_low = high_verts_of_lows[j];
    unsigned high_buf[MAX_UP];
    unsigned high_buf_size = 0;
    for (unsigned k = 0; k < verts_per_low; ++k) {
      unsigned vert = verts_of_high[high_verts_of_low[k]];
      unsigned first_use = lows_of_verts_offsets[vert];
      unsigned end_use = lows_of_verts_offsets[vert + 1];
      if (k) {
        high_buf_size = intersect(
            high_buf,
            high_buf_size,
            lows_of_verts + first_use,
            end_use - first_use);
      } else if (dual_mode) {
        assert(end_use - first_use <= MAX_UP);
        high_buf_size = copy_except(
            lows_of_verts + first_use,
            high_buf,
            end_use - first_use,
            i);
      } else {
        assert(end_use - first_use <= MAX_UP);
        high_buf_size = copy(
            lows_of_verts + first_use,
            high_buf,
            end_use - first_use);
      }
    }
    if (dual_mode) {
      assert(high_buf_size <= 1);
      lows_of_high[j] = ((high_buf_size) ? high_buf[0] : INVALID);
    } else {
      assert(high_buf_size == 1);
      lows_of_high[j] = high_buf[0];
    }
  }
}

/* TODO: This is the #1 most expensive function, takes up 70% of
   adaptation time !
   If you are going to optimize, optimize here !
   1) try to make this marginally faster
   2) look for higher-level ways to avoid calling this function
 */
static unsigned* reflect_down_general(
    unsigned dual_mode,
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* verts_of_highs,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts)
{
  unsigned verts_per_high = the_down_degrees[high_dim][0];
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned verts_per_low = the_down_degrees[low_dim][0];
  unsigned* lows_of_highs = LOOP_MALLOC(unsigned, nhighs * lows_per_high);
  unsigned const* const* high_verts_of_lows =
    the_canonical_orders[high_dim][low_dim][0];
  LOOP_EXEC(reflect_down_entity, nhighs,
      dual_mode,
      verts_of_highs,
      lows_of_verts_offsets,
      lows_of_verts,
      verts_per_high,
      lows_per_high,
      verts_per_low,
      lows_of_highs,
      high_verts_of_lows);
  return lows_of_highs;
}

/* some optimization here due to the relatively
   high cost of these operations.
   first, nverts_wanted is either 2 or 3.
   that is because reflect_down is only
   used to derive intermediate downward adjacencies,
   and we always have the n->0 adjacency, so it must
   be one of:
     3->2
     3->1
     2->1
   and for d=2, d=1, the number of vertices of a simplex
   is n=3, n=2, respectively.
   we can use this macro (template) to make it clearer to
   the compiler it should generate unrolled code for
   these two cases.
*/

#define FIND_LOW_FAST(N) \
LOOP_INOUT static inline unsigned find_low_fast_##N( \
    unsigned const* verts_wanted, \
    unsigned const* verts_of_lows, \
    unsigned const* lows_of_verts_offsets, \
    unsigned const* lows_of_verts) \
{ \
  unsigned v0 = verts_wanted[0]; \
  unsigned f = lows_of_verts_offsets[v0]; \
  unsigned e = lows_of_verts_offsets[v0 + 1]; \
  for (unsigned i = f; i < e; ++i) { \
    unsigned low = lows_of_verts[i]; \
    unsigned const* verts_of_low = verts_of_lows + low * N; \
    unsigned j; \
    for (j = 0; j < N; ++j) { \
      unsigned wanted_vert = verts_wanted[j]; \
      unsigned k; \
      for (k = 0; k < N; ++k) \
        if (verts_of_low[k] == wanted_vert) \
          break; \
      if (k == N) \
        break; \
    } \
    if (j == N) \
      return low; \
  } \
  assert(0); \
}

FIND_LOW_FAST(2)
FIND_LOW_FAST(3)

LOOP_INOUT static inline unsigned find_low_fast(
    unsigned nverts_wanted,
    unsigned const* verts_wanted,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts)
{
  if (nverts_wanted == 2)
    return find_low_fast_2(verts_wanted, verts_of_lows, lows_of_verts_offsets,
        lows_of_verts);
  else
    return find_low_fast_3(verts_wanted, verts_of_lows, lows_of_verts_offsets,
        lows_of_verts);
}

LOOP_KERNEL(reflect_down_entity_fast,
    unsigned const* verts_of_highs,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts,
    unsigned verts_per_high,
    unsigned lows_per_high,
    unsigned verts_per_low,
    unsigned* lows_of_highs,
    unsigned const* const* high_verts_of_lows)

  unsigned const* verts_of_high = verts_of_highs + i * verts_per_high;
  unsigned* lows_of_high = lows_of_highs + i * lows_per_high;
  for (unsigned j = 0; j < lows_per_high; ++j) {
    unsigned const* high_verts_of_low = high_verts_of_lows[j];
    unsigned verts_wanted[3];
    for (unsigned k = 0; k < verts_per_low; ++k)
      verts_wanted[k] = verts_of_high[high_verts_of_low[k]];
    lows_of_high[j] = find_low_fast(verts_per_low, verts_wanted,
        verts_of_lows, lows_of_verts_offsets, lows_of_verts);
  }
}

static unsigned* reflect_down_fast(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* verts_of_highs,
    unsigned const* verts_of_lows,
    unsigned const* lows_of_verts_offsets,
    unsigned const* lows_of_verts)
{
  unsigned verts_per_high = the_down_degrees[high_dim][0];
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned verts_per_low = the_down_degrees[low_dim][0];
  unsigned* lows_of_highs = LOOP_MALLOC(unsigned, nhighs * lows_per_high);
  unsigned const* const* high_verts_of_lows =
    the_canonical_orders[high_dim][low_dim][0];
  LOOP_EXEC(reflect_down_entity_fast, nhighs,
      verts_of_highs,
      verts_of_lows,
      lows_of_verts_offsets,
      lows_of_verts,
      verts_per_high,
      lows_per_high,
      verts_per_low,
      lows_of_highs,
      high_verts_of_lows);
  return lows_of_highs;
}

unsigned* mesh_reflect_down(
    struct mesh* m,
    unsigned high_dim,
    unsigned low_dim)
{
  unsigned nhighs = mesh_count(m, high_dim);
  unsigned const* verts_of_highs = mesh_ask_down(m, high_dim, 0);
  struct const_up* lows_of_verts = mesh_ask_up(m, 0, low_dim);
  unsigned const* verts_of_lows = mesh_ask_down(m, low_dim, 0);
  {
    double lows_per_high = mesh_estimate_degree(m, high_dim, low_dim);
    double verts_per_low = mesh_estimate_degree(m, low_dim, 0);
    double lows_per_vert = mesh_estimate_degree(m, 0, low_dim);
    printf("mesh_reflect_down(high=%u, low=%u)\n", high_dim, low_dim);
    double sets_per_high = lows_per_high * verts_per_low;
    printf("set queries/intersections per high = (%f*%f) = %f\n",
        lows_per_high, verts_per_low, sets_per_high);
    double set_size = lows_per_vert;
    printf("average set size %f\n", set_size);
    printf("best case per high (intersect cost O(N)) = %f\n",
        sets_per_high * set_size);
    printf("worst case per high (intersect cost O(N^2)) = %f\n",
        sets_per_high * set_size * set_size);
  }
  {
    printf("predicted improvement:\n");
    double lows_per_high = mesh_estimate_degree(m, high_dim, low_dim);
    double verts_per_low = mesh_estimate_degree(m, low_dim, 0);
    double lows_per_vert = mesh_estimate_degree(m, 0, low_dim);
    double num_sets = lows_per_high * lows_per_vert;
    printf("number of sets %f\n", num_sets);
    double set_size = verts_per_low;
    printf("average set size %f\n", set_size);
    printf("cost per high: %f\n", num_sets * set_size * set_size);
  }
  double t0 = MPI_Wtime();
//unsigned* out = reflect_down(high_dim, low_dim, nhighs,
//    verts_of_highs, lows_of_verts->offsets, lows_of_verts->adj);
  unsigned* out = reflect_down_fast(high_dim, low_dim, nhighs,
      verts_of_highs, verts_of_lows, lows_of_verts->offsets, lows_of_verts->adj);
  double t1 = MPI_Wtime();
  printf("actual cost: %f (total) %e (per high)\n", t1 - t0,
      (t1 - t0)/((double)nhighs));
  return out;
}

unsigned* get_dual(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts)
{
  return reflect_down_general(1, elem_dim, elem_dim - 1, nelems, verts_of_elems,
      elems_of_verts_offsets, elems_of_verts);
}

unsigned* mesh_get_dual(struct mesh* m)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  struct const_up* elems_of_verts = mesh_ask_up(m, 0, elem_dim);
  return get_dual(elem_dim, nelems, verts_of_elems,
      elems_of_verts->offsets, elems_of_verts->adj);
}
