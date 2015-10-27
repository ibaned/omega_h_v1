#include "invert_map.h"

#include "ints.h"
#include "loop.h"


__global__ Worker_Kernal( unsigned * in ,
		unsigned * offsets , unsigned * counts,
		unsigned *out )
{
	int i = (blockIdx.x * blockDim.x + threadIdx.x);
	out[ offsets[in[i]] + counts[in[i]]++  ] = i; //I mean....
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
    counts[in[i]]++; //Broken Operation.
  unsigned* offsets = uints_exscan(counts, nout);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  uints_zero(counts, nout);
  for (unsigned i = 0; i < nin; ++i) {
    unsigned d = in[i];
    unsigned o = offsets[d];
    unsigned j = counts[d]++;
    out[o + j] = i;
  }
  loop_free(counts);
  *p_out = out;
  *p_offsets = offsets;
}
