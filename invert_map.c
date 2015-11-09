#include "invert_map.h"

#include "ints.h"
#include "loop.h"

LOOP_KERNEL(fill, unsigned const* in,
		unsigned* offsets, unsigned* counts,
		unsigned* out)
  unsigned d = in[i];
  unsigned o = offsets[d];
  unsigned j = loop_atomic_increment(&(counts[d]));
  out[o + j] = i;
}

/* the atomic increments give the fill
   kernel some non-determinism,
   which we can counteract by sorting the sub-arrays.
   given the assumption that these sub-arrays are small
   (about a dozen entries),
   we use a hand-coded selection sort.

   there is a reward for someone who comes up with
   an equally efficient implementation that is deterministic
   from the start */
LOOP_KERNEL(sort,
		unsigned* offsets,
		unsigned* out)
  unsigned first = offsets[i];
  unsigned end = offsets[i + 1];
  for (unsigned j = first; j < end; ++j) {
    unsigned min_k = j;
    for (unsigned k = j + 1; j < end; ++j)
      if (out[k] < out[min_k])
        min_k = k;
    unsigned tmp = out[j];
    out[j] = out[min_k];
    out[min_k] = tmp;
  }
}

void invert_map(
    unsigned nin,
    unsigned const* in,
    unsigned nout,
    unsigned** p_out,
    unsigned** p_offsets)
{
  unsigned* counts = LOOP_MALLOC(unsigned, nout);
  uints_zero(counts, nout);
  for (unsigned i = 0; i < nin; ++i)
    counts[in[i]]++; /* TODO: loop kernel for this */
  unsigned* offsets = uints_exscan(counts, nout);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  uints_zero(counts, nout);
  LOOP_EXEC(fill, nin, in, offsets, counts, out);
  loop_free(counts);
  LOOP_EXEC(sort, nout, offsets, out);
  *p_out = out;
  *p_offsets = offsets;
}
