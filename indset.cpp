#include "indset.hpp"

#include <cstdlib>

#include "arrays.hpp"
#include "comm.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"

/* given an undirected graph (usually obtained from get_star),
 * an initial filtered subset of the vertices,
 * and a scalar value (goodness) at the vertices,
 * returns a maximal independent set of the subgraph
 * induced by the filter, whose members are
 * preferrably those with locally high goodness values.
 * ties in goodness values are broken by choosing the
 * vertex with the lowest "global" value.
 * In order for this function to run efficiently,
 * there should not exist long paths with monotonically
 * decreasing goodness.
 */

/* the runtime of the independent set algorithm
 * as written below is O(iterations * vertices).
 * as such we would want iterations to be small
 * and this constant is a hard ceiling above
 * which the function aborts the program
 */
#define MAX_ITERATIONS 100

enum {
  NOT_IN_SET = 0,
  IN_SET = 1,
  UNKNOWN = 2
};

/* this is a modified version of Luby's Maximal Independent Set
   algorithm, published in:

   Luby, Michael.
   "A simple parallel algorithm for the maximal independent set problem."
   SIAM journal on computing 15.4 (1986): 1036-1053.

   instead of random values at the vertices, we use "goodness" values.
   these values are usually the quality of mesh cavities, which
   have the similar desirable property that there are no
   long paths in the graph with monotonically increasing goodness.
   (the algorithm's iteration count is proportional to the length
    of the longest such path). */

LOOP_KERNEL(indset_at_vert,
    unsigned const* offsets,
    unsigned const* adj,
    double const* goodness,
    unsigned long const* global,
    unsigned const* old_state,
    unsigned* state)

  if (old_state[i] != UNKNOWN)
    return;
  unsigned first_adj = offsets[i];
  unsigned end_adj = offsets[i + 1];
  for (unsigned j = first_adj; j < end_adj; ++j)
    if (old_state[adj[j]] == IN_SET) {
      state[i] = NOT_IN_SET;
      return;
    }
  double myg = goodness[i];
  for (unsigned j = first_adj; j < end_adj; ++j) {
    unsigned other = adj[j];
    if (old_state[other] == NOT_IN_SET)
      continue;
    double og = goodness[other];
    if (myg == og && global[other] < global[i])
      return;
    if (myg < og)
      return;
  }
  state[i] = IN_SET;
}

LOOP_KERNEL(init_state,
    unsigned const* filter,
    unsigned* state)
  if (filter[i])
    state[i] = UNKNOWN;
  else
    state[i] = NOT_IN_SET;
}

static unsigned* find_indset(
    struct mesh* m,
    unsigned ent_dim,
    unsigned nverts,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned const* filter,
    double const* goodness,
    unsigned long const* global)
{
  unsigned* state = LOOP_MALLOC(unsigned, nverts);
  LOOP_EXEC(init_state, nverts, filter, state);
  for (unsigned it = 0; it < MAX_ITERATIONS; ++it) {
    unsigned* old_state = uints_copy(state, nverts);
    LOOP_EXEC(indset_at_vert, nverts,
        offsets, adj, goodness, global, old_state, state);
    loop_free(old_state);
    mesh_conform_uints(m, ent_dim, 1, &state);
    if (comm_max_uint(uints_max(state, nverts)) < UNKNOWN)
      return state;
  }
  LOOP_NORETURN(0);
}

unsigned* mesh_find_indset(struct mesh* m, unsigned ent_dim,
    unsigned const* candidates, double const* qualities)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nents = mesh_count(m, ent_dim);
  unsigned const* star_offsets =
    mesh_ask_star(m, ent_dim, elem_dim)->offsets;
  unsigned const* star =
    mesh_ask_star(m, ent_dim, elem_dim)->adj;
  unsigned long const* global = 0;
  unsigned long* to_free = 0;
  if (mesh_is_parallel(m)) {
    global = mesh_ask_globals(m, ent_dim);
  } else {
    global = to_free = ulongs_linear(nents, 1);
  }
  unsigned* indset = find_indset(m, ent_dim, nents,
      star_offsets, star, candidates, qualities, global);
  loop_free(to_free);
  return indset;
}

unsigned* mesh_indset_offsets(struct mesh* m, unsigned ent_dim,
    unsigned const* candidates, double const* qualities)
{
  unsigned nents = mesh_count(m, ent_dim);
  unsigned* indset = mesh_find_indset(m, ent_dim, candidates, qualities);
  unsigned* offsets = uints_exscan(indset, nents);
  loop_free(indset);
  return offsets;
}
