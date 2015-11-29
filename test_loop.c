#include "loop.h"
#include "ints.h"
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <stdlib.h>

#define N 100000000

static void Ints_Test()
{
  unsigned* a;
  unsigned* b;
  unsigned* c;
  unsigned* d;
  unsigned zz;
  clock_t start;
  clock_t end;

  //Zero (Runtime)
  a = uints_filled(N, 0);
  start = clock();
  end = clock();
  b = (unsigned*) loop_to_host(a, N * sizeof(unsigned));
  for(unsigned i = 0; i < N; ++i)
    if(b[i] != 0) {
      printf("zero failed");
      abort();
    }
  loop_host_free(b);
  loop_free(a);
  printf("ZERO: %F\n", (double)(end-start)/CLOCKS_PER_SEC);

  //units_copy (Max values ?? ; Min values ?? ; Runtime)
  a = LOOP_HOST_MALLOC(unsigned, N);
  for(unsigned i = 0; i < N ;i++)
  {
    a[i] = ((unsigned)rand()) % 0xFFFFFFFF;
  }
  c = (unsigned*) loop_to_device(a, N * sizeof(unsigned));
  start = clock();
  d = uints_copy(c, N);
  end = clock();
  b = (unsigned*) loop_to_host(d, N * sizeof(unsigned));
  loop_free(c);
  loop_free(d);
  for(unsigned i = 0; i < N; ++i)
    if(a[i] != b[i]) {
      printf("uints copy failed!");
      abort();
    }
  loop_host_free(b);
  loop_host_free(a);
  printf( "COPY: %F\n" , (double)(end-start)/CLOCKS_PER_SEC);

  //units_max (Max values; Min values; Runtime)
  a = (unsigned*) malloc(N * sizeof(unsigned));
  for(unsigned i = 0; i < N; ++i)
  {
    a[i] = ((unsigned)rand()) % 0xFFFFFFFF;
  }
  a[0] = 0xFFFFFFFF; //Max value
  c = (unsigned*) loop_to_device(a, N * sizeof(unsigned));
  start = clock();
  zz = uints_max(c, N);
  end = clock();
  loop_free(c);
  loop_host_free(a);
  if (zz != 0xFFFFFFFF) {
    printf( "MAX FAILED");
    abort();
  }
  printf( "MAX: %F\n" , (double)(end-start)/CLOCKS_PER_SEC);

  // uints_exscan (Max values; Min values; Runtime)
  a = LOOP_HOST_MALLOC(unsigned, N);
  for( int i =0 ; i< N ;i++)
  {
    a[i] = ((unsigned)rand()) % 0xFFFFFFFF;
  }
  /* todo ? */
  loop_host_free(a);

  // uints_unscan (Max values; Min values; Runtime)

  // uints_negate (I'm confused)
}

LOOP_KERNEL(count,
    unsigned const* in,
    unsigned* counts)
  loop_atomic_increment(&(counts[in[i]]));
}

static void COUNT_TEST()
{
  unsigned  a[7] = {1,2,1,3,5,1,2};
  unsigned* b = (unsigned*) loop_to_device(a, sizeof(unsigned) * 7);

  unsigned* counts = uints_filled(7, 0);

  LOOP_EXEC(count, 7, b, counts);
  loop_free(b);

  unsigned* z = (unsigned *)loop_to_host(counts, sizeof(unsigned) * 7);

  for(int i = 0 ; i<7; i++)
  {
    printf("%d\n", z[i]);
  }
  loop_free(counts);
  loop_host_free(z);
}


int main()
{
#ifdef __CUDACC__
  cudaDeviceReset();
#endif

  unsigned a[7] = {1,2,1,3,5,1,2};
  unsigned* b = (unsigned*) loop_to_device(a, sizeof(unsigned) * 7);

  unsigned* counts = uints_filled(7,0);

  LOOP_EXEC(count, 7, b, counts);
  loop_free(b);

  unsigned* z = (unsigned*) loop_to_host(counts, sizeof(unsigned) * 7);

  for(unsigned i = 0; i < 7; ++i)
  {
    printf("%d\n", z[i]);
  }
  loop_free(counts);
  loop_host_free(z);

  Ints_Test();
  COUNT_TEST();
  return 0;
}
