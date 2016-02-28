#include "invert_map.h"


#include "arrays.h"
#include "ints.h"
#include "loop.h"

#include <thrust/sort.h>

struct Counter
{
  unsigned origin;
  unsigned count;
};

LOOP_KERNEL(count,
    unsigned const* in,
    unsigned* counts)
  loop_atomic_increment(&(counts[in[i]]));
}

LOOP_KERNEL(fill, unsigned const* in,
    unsigned* offsets,
    unsigned* counts,
    unsigned* out,
    Counter* counter)
  unsigned d = in[i];
  unsigned o = offsets[d];
  unsigned j = counts[d]+ counter[i].count;
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




LOOP_KERNEL( assign , struct Counter * out )
  out[i].origin = i;
}

LOOP_KERNEL( Pointer_assign , struct Counter * out, struct Counter** ref )
  ref[i] = &(out[i]);
}

LOOP_KERNEL( count_work , struct Counter ** ref , unsigned * in )
  ref[i]->count = 0;
  __syncthreads();
  while( i > 0 && in[i]==in[i-1])
  {
	atomicAdd( &(ref[i]->count) , 1 );
	i--;
  }
}

LOOP_KERNEL( copy , unsigned const* in , unsigned * out)
  out[i] = in[i];
}


#include <thrust/device_vector.h>
#include <thrust/fill.h>

void Count_Sort_Dance(unsigned const*in , struct Counter*  out, unsigned nin)
{
  struct Counter** ref = LOOP_MALLOC( struct Counter* , nin);
  LOOP_EXEC( assign, nin , out );
  unsigned * n_sorted = LOOP_MALLOC( unsigned , nin);
  //unsigned *a = (unsigned*)malloc( sizeof(int) * nin);
  LOOP_EXEC( copy , nin , in , n_sorted);
  LOOP_EXEC( Pointer_assign, nin, out, ref);
  //LOOP_EXEC( Pointer_assign, nin, out, ref);
  //thrust::device_ptr<unsigned> t_n_sorted(n_sorted);
  //thrust::device_ptr<unsigned> t_n_sorted_end(n_sorted+nin);
  //thrust::device_ptr<struct Counter*> t_ref(ref);
  //thrust::device_vector<unsigned>

  //thrust::device_vector<unsigned> a_sorted( n_sorted, n_sorted+nin);
  //thrust::device_vector<Counter*>
  /*
  thrust::fill(thrust::device_ptr<unsigned> (n_sorted),
		  thrust::device_ptr<unsigned> (n_sorted+nin),
		  0);
  printf(" DID I DIE!?!?!?\n");
  */

  cudaError_t ret = cudaGetLastError();
  assert(ret == cudaSuccess);

  thrust::stable_sort_by_key(
    thrust::device_ptr<unsigned> (n_sorted) ,
    thrust::device_ptr<unsigned> (n_sorted+nin),
    thrust::device_ptr<struct Counter*> (ref) );

  //printf("I'M NOT DEAD\n");

  LOOP_EXEC( count_work , nin, ref , n_sorted);
  loop_free(ref);

  loop_free(n_sorted);
}




void invert_map(
    unsigned nin,
    unsigned const* in,
    unsigned nout,
    unsigned** p_out,
    unsigned** p_offsets)
{
  unsigned* counts = uints_filled(nout, 0);
  LOOP_EXEC(count, nin, in, counts);
  struct Counter* counters = LOOP_MALLOC( struct Counter , nin);
  Count_Sort_Dance( in, counters , nin);
  unsigned* offsets = uints_exscan(counts, nout);
  unsigned* out = LOOP_MALLOC(unsigned, nin);
  loop_free(counts);
  counts = uints_filled(nout, 0);
  LOOP_EXEC(fill, nin, in, offsets, counts, out, counters);
  loop_free(counts);
  //printf("BEFORE:\n");
  //int i = 0;
  //for(i =0 ; i< nout ; i++) printf("%i : %i \n", i , out[i]);
  LOOP_EXEC(sort, nout, offsets, out);
  //printf("AFTER:\n");
  //for(i =0 ; i< nout ; i++) printf("%i : %i \n", i , out[i]);
  loop_free(counters);
  //printf("DId this kill me\n");
  *p_out = out;
  *p_offsets = offsets;
}
