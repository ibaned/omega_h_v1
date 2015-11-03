#include "invert_map.h"

#include "ints.h"
#include "loop.h"

LOOP_KERNEL(fill, unsigned const* in,
	unsigned* offsets,
	unsigned* counts,
	unsigned* out)
  unsigned d = in[i];
  unsigned o = offsets[d];
  unsigned j = loop_atomic_increment(&(counts[d]));
  out[o + j] = i;
}

LOOP_KERNEL(count ,
	unsigned const* in ,
	unsigned * counts)
  loop_atomic_increment(&(counts[in[i]]));
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
  LOOP_EXEC(count, nin, in, counts);
  unsigned* offsets = uints_exscan(counts, nout);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  uints_zero(counts, nout);
  LOOP_EXEC(fill, nin, in, offsets, counts, out);
  loop_free(counts);
  *p_out = out;
  *p_offsets = offsets;
}
