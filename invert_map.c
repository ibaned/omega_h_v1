#include "invert_map.h"

#include "ints.h"
#include "loop.h"

#include <thrust/sort.h>


LOOP_KERNEL(fill, unsigned const* in,
  unsigned* offsets,
  unsigned* counts,
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
   we use a hand-coded selection sort.git
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
    for (unsigned k = j + 1; k < end; ++k)
      if (out[k] < out[min_k])
        min_k = k;
    unsigned tmp = out[j];
    out[j] = out[min_k];
    out[min_k] = tmp;
  }
}

LOOP_KERNEL(count ,
  unsigned const* in ,
  unsigned * counts)
  loop_atomic_increment(&(counts[in[i]]));
}


struct Counter
{
  unsigned origin;
  unsigned sorted;
  unsigned count;
};

LOOP_KERNEL( assign , struct Counter * out )
  out[i].origin = i;
}

LOOP_KERNEL( count_work , struct Counter * out , unsigned * in )
  out[i].count = 0;
  while( i > 0 && in[i]==in[i-1])
  {
	  out[i].count += 1;
	  i--;
  }
}

void Count_Sort_Dance( unsigned *in , struct Counter*  out, unsigned nin)
{
  out = LOOP_MALLOC( struct Counter , nin);
  LOOP_EXEC( assign, nin , out );
  unsigned * n_sorted = uints_copy( in , nin);

  thrust::stable_sort_by_key( n_sorted ,n_sorted+nin , out );
  LOOP_EXEC( assign, nin, out);

  LOOP_EXEC( count_work , nin, out , in);
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
  LOOP_EXEC(sort, nout, offsets, out);
  *p_out = out;
  *p_offsets = offsets;
}
