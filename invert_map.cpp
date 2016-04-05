#include "invert_map.hpp"

#include <cassert>
#include <chrono>
#include <algorithm>

#include "arrays.hpp"
#include "ints.hpp"
#include "loop.hpp"

#ifdef LOOP_CUDA_HPP
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#endif

namespace omega_h {

unsigned long invert_map_total_nin = 0;
double invert_map_time = 0;
double invert_map_count_time = 0;
double invert_map_fill_time = 0;
double invert_map_sort_time = 0;

#if 0

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
  auto t0 = std::chrono::high_resolution_clock::now();
  struct Counter* counters = LOOP_MALLOC(Counter, nin);
  LOOP_EXEC(assign, nin, in, counters);
#ifdef LOOP_CUDA_HPP
  thrust::stable_sort(thrust::device_ptr<Counter>(counters),
      thrust::device_ptr<Counter>(counters + nin));
#else
  std::sort(counters, counters + nin);
#endif
  unsigned* aoffsets = filled_array<unsigned>(nout + 1, 0);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  LOOP_EXEC(fill,nin, aoffsets, out, counters);
  array_set(aoffsets, nout, nin);
  loop_free(counters);
  *p_out = out;
  *p_offsets = aoffsets;
  auto t1 = std::chrono::high_resolution_clock::now();
  auto diff = t1 - t0;
  invert_map_time += double(std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count() * 1e-9);
  invert_map_total_nin += nin;
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

#if 0
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
#endif

void invert_map(
    unsigned nin,
    unsigned const* in,
    unsigned nout,
    unsigned** p_out,
    unsigned** p_offsets)
{
  auto t0 = std::chrono::high_resolution_clock::now();
  unsigned* counts = filled_array<unsigned>(nout, 0);
#ifdef LOOP_CUDA_HPP
  CUDACALL(cudaDeviceSynchronize());
#endif
  auto c0 = std::chrono::high_resolution_clock::now();
  LOOP_EXEC(count, nin, nout, in, counts);
#ifdef LOOP_CUDA_HPP
  CUDACALL(cudaDeviceSynchronize());
#endif
  auto c1 = std::chrono::high_resolution_clock::now();
  auto diff = c1 - c0;
  invert_map_count_time += double(std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count() * 1e-9);
  unsigned* offsets = uints_exscan(counts, nout);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  loop_free(counts);
  counts = filled_array<unsigned>(nout, 0);
  auto f0 = std::chrono::high_resolution_clock::now();
  LOOP_EXEC(fill, nin, in, offsets, counts, out);
#ifdef LOOP_CUDA_HPP
  CUDACALL(cudaDeviceSynchronize());
#endif
  auto f1 = std::chrono::high_resolution_clock::now();
  diff = f1 - f0;
  invert_map_fill_time += double(std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count() * 1e-9);
  loop_free(counts);
#if 0
  auto s0 = std::chrono::high_resolution_clock::now();
  LOOP_EXEC(sort, nout, offsets, out);
#ifdef LOOP_CUDA_HPP
  CUDACALL(cudaDeviceSynchronize());
#endif
  auto s1 = std::chrono::high_resolution_clock::now();
  diff = s1 - s0;
  invert_map_sort_time += double(std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count() * 1e-9);
#endif
  *p_out = out;
  *p_offsets = offsets;
  auto t1 = std::chrono::high_resolution_clock::now();
  diff = t1 - t0;
  invert_map_time += double(std::chrono::duration_cast<std::chrono::nanoseconds>(diff).count() * 1e-9);
  invert_map_total_nin += nin;
}

#endif

}
