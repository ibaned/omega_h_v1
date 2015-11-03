#include "loop.h"
#include "ints.h"
#include <stdio.h>
LOOP_KERNEL(fill, double* a)
  a[i] = 0;
}

LOOP_KERNEL(count ,
	unsigned const* in ,
	unsigned * counts)
	loop_atomic_increment(&(counts[in[i]]));
}

int main()
{

	/*
  unsigned n = 1024*1024;
  double* d = LOOP_MALLOC(double, n);
  LOOP_EXEC(fill, n, d);
  loop_free(d);
*/
  cudaDeviceReset();

  unsigned  a[7] = {1,2,1,3,5,1,2};
  unsigned *b  = (unsigned *)loop_to_device(a , sizeof(unsigned)*7 );

  unsigned* counts = LOOP_MALLOC(unsigned, 7);
  uints_zero(counts, 7);

  LOOP_EXEC( count , 7 , b , counts);

  unsigned *z = (unsigned *)loop_cuda_to_host( counts, sizeof(unsigned)*7);

  for(int i = 0 ; i<7; i++)
  {
	  printf( "%d\n" ,z[i]);
  }
  //loop_free( counts);

}
