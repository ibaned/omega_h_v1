#include "invert_map.h"

#include "arrays.h"
#include "ints.h"
#include "loop.h"

#ifdef LOOP_CUDA_H

#include <thrust/sort.h>
#include <thrust/device_ptr.h>

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

LOOP_KERNEL(assign, const unsigned* in, Counter* ref )
  ref[i].origin = i;
  ref[i].value  = in[i];
}

LOOP_IN bool operator<(const Counter &A, const Counter &B)
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
  struct Counter* counters = LOOP_MALLOC(Counter, nin);
  LOOP_EXEC(assign, nin , in, counters );
  thrust::sort(
      thrust::device_ptr<Counter>(counters),
      thrust::device_ptr<Counter>(counters + nin));
  unsigned* aoffsets =uints_filled(nout + 1, 0);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  LOOP_EXEC(fill, nin,aoffsets, out, counters);
  CUDACALL(cudaMemcpy(aoffsets + nout, &(nin),
        sizeof(unsigned), cudaMemcpyHostToDevice));
  loop_free(counters);
  *p_out = out;
  *p_offsets = aoffsets;
}

#else

LOOP_KERNEL(count,
    unsigned nout,
    unsigned const* in,
    unsigned* counts)
  unsigned d = in[i];
  assert(d < nout);
  loop_atomic_increment(&(counts[d]));
}

LOOP_KERNEL(fill,
    unsigned const* in,
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
    for (unsigned k = j + 1; k < end; ++k)
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
  unsigned* counts = uints_filled(nout, 0);
  LOOP_EXEC(count, nin, nout, in, counts);
  unsigned* offsets = uints_exscan(counts, nout);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  loop_free(counts);
  counts = uints_filled(nout, 0);
  LOOP_EXEC(fill, nin, in, offsets, counts, out);
  loop_free(counts);
  LOOP_EXEC(sort, nout, offsets, out);
  *p_out = out;
  *p_offsets = offsets;
}

#endif
