#include "invert_map.h"

#include "arrays.h"
#include "ints.h"
#include "loop.h"

#include <thrust/sort.h>

struct Counter
{
  unsigned origin;
  unsigned value;
};

LOOP_KERNEL(fill,
    unsigned* offset,
    unsigned* out,
    Counter* counter)
  out[i] = counter[i].origin;
  if(i!=0)
  {
    if( counter[i].value != counter[i-1].value)
    {
      offset[counter[i].value] = i;
    }
  }
}

LOOP_KERNEL( assign ,const unsigned * in, struct Counter * ref )
  ref[i].origin = i;
  ref[i].value  = in[i];
}

#include <thrust/device_vector.h>
#include <thrust/fill.h>

__host__ __device__ bool operator<(const Counter &A, const Counter &B)
{
  return (A.value < B.value);
}

void invert_map(
    unsigned nin,
    unsigned const* in,
    unsigned nout,
    unsigned** p_out,
    unsigned** p_offsets)
{
  struct Counter* counters = LOOP_MALLOC( struct Counter, nin);
  LOOP_EXEC( assign, nin , in, counters );

  thrust::sort(
      thrust::device_ptr<struct Counter> (counters),
      thrust::device_ptr<struct Counter> (counters+nin));

  unsigned* aoffsets =uints_filled(nout+1, 0);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  unsigned i = nin;
  LOOP_EXEC(fill, nin,aoffsets, out, counters);
  CUDACALL(cudaMemcpy(aoffsets + nout, &(i), sizeof(unsigned), cudaMemcpyHostToDevice));
  loop_free(counters);
  *p_out = out;
  *p_offsets = aoffsets;
}
