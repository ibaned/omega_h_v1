#include "indset.h"

#include <stdlib.h>

#include "ints.h"
#include "loop.h"

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

static void at_vert(
    unsigned const* offsets,
    unsigned const* adj,
    double const* goodness,
    unsigned const* old_state,
    unsigned* state,
    unsigned i)
{
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
    if (old_state[adj[j]] == NOT_IN_SET)
      continue;
    double og = goodness[adj[j]];
    if (myg == og && adj[j] < i)
      return;
    if (myg < og)
      return;
  }
  state[i] = IN_SET;
}

unsigned* find_indset(
    unsigned nverts,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned const* filter,
    double const* goodness)
{
  unsigned* state = LOOP_MALLOC(unsigned, nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    if (filter[i])
      state[i] = UNKNOWN;
    else
      state[i] = NOT_IN_SET;
  }
  for (unsigned it = 0; it < MAX_ITERATIONS; ++it) {
    unsigned* old_state = uints_copy(state, nverts);
    for (unsigned i = 0; i < nverts; ++i)
      at_vert(offsets, adj, goodness, old_state, state, i);
    LOOP_FREE(old_state);
    if (uints_max(state, nverts) < UNKNOWN)
      return state;
  }
  abort();
}
