#include "loop.h"
#include "ints.h"
#include <stdio.h>
#include <time.h>
#include <limits.h>

#define N 1000000




LOOP_KERNEL(fill, double* a)
  a[i] = 0;
}

LOOP_KERNEL(count ,
	unsigned const* in ,
	unsigned * counts)
	loop_atomic_increment(&(counts[in[i]]));
}

void Ints_Test()
{
  unsigned *a, *b, *c , *d , zz; clock_t start , end;
  //Zero (Runtime)
  a = LOOP_MALLOC(unsigned , N);
  start = clock();
  uints_zero( a , N);
  end = clock();

#ifdef __CUDACC__
  b = (unsigned *)loop_cuda_to_host( a , N * sizeof(unsigned));
  for( unsigned i = 0 ; i < N ; i++)
    if( b[i] != 0){ printf("zero failed"); abort();}
  free(b);
#else
  for( unsigned i = 0 ; i < N ; i++)
	  if( a[i] != 0){ printf("zero failed"); abort();}
#endif

  loop_free(a);;
  printf( "ZERO: %f" , (double)(end-start)/CLOCKS_PER_SEC);


  //units_copy (Max values ?? ; Min values ?? ; Runtime)
  a = (unsigned *)malloc(N *sizeof(unsigned));
  for( int i =0 ; i< N ;i++)
  {
	  a[i] = rand()%(0xFFFFFFFF);
  }
#ifdef __CUDACC__
  c = (unsigned *)loop_cuda_to_device( a , N* sizeof(unsigned));
  start = clock();
  d = uints_copy ( c , N);
  end = clock();
  b = (unsigned *)loop_cuda_to_host( d , N*sizeof(unsigned));

  loop_free(c); loop_free(d);
#else
  start = clock();
  b = uints_copy( a , N);
  end = clock();
#endif

  for(unsigned i = 0; i<N;i++)
    if( a[i] != b[i]){ printf("uints copy failed!"); abort();}
  free(b); free(a);
  printf( "COPY: %f" , (double)(end-start)/CLOCKS_PER_SEC);


  //units_max (Max values; Min values; Runtime)

  a = (unsigned *)malloc(N *sizeof(unsigned));
  for( int i =0 ; i< N ;i++)
  {
	  a[i] = rand()%(0xFFFFFFFF);
  }
  a[0] = 0xFFFFFFFF; //Max value

#ifdef __CUDACC__
  c = (unsigned *)loop_cuda_to_device( a , N* sizeof(unsigned));
  start = clock();
  zz = uints_max ( c , N);
  end = clock();
  loop_free(c);
#else
  start = clock();
  zz = uints_max( a , N);
  end = clock();
#endif
  free(a);

  if ( zz != 0xFFFFFFFF){ printf( "MAX FAILED"); abort();}

  printf( "COPY: %f" , (double)(end-start)/CLOCKS_PER_SEC);

  // uints_exscan (Max values; Min values; Runtime)


  // uints_unscan (Max values; Min values; Runtime)




  // uints_negate (I'm confused)



}

int COUNT_TEST()
{
	unsigned  a[7] = {1,2,1,3,5,1,2};
	  unsigned *b = (unsigned *)loop_to_device(a , sizeof(unsigned)*7 );

	  unsigned* counts = LOOP_MALLOC(unsigned, 7);
	  uints_zero(counts, 7);

	  LOOP_EXEC( count , 7 , b , counts);

	  unsigned *z = (unsigned *)loop_cuda_to_host( counts, sizeof(unsigned)*7);

	  for(int i = 0 ; i<7; i++)
	  {
		  printf( "%d\n" ,z[i]);
	  }
	  loop_free( counts);

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
  loop_free( counts);

}
