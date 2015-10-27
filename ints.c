#include "ints.h"

#include "loop.h"

#ifdef __CUDACC__
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include "mycuda.cuh"
#include <thrust/functional.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#endif



#ifdef __CUDACC__
void uints_zero(unsigned* a, unsigned n)
{
	thrust::device_ptr< unsigned int> p (a);
	thrust::fill( p , p+n , (unsigned) 0);

}

unsigned* uints_copy(unsigned  * a , unsigned n)
{
  unsigned *b = LOOP_MALLOC(unsigned , n );
  CUDACALL( cudaMemcpy(b , a , n*sizeof(unsigned) , cudaMemcpyDeviceToDevice));

  return b;
}

unsigned* uints_copy(unsigned const * a , unsigned n)
{
  unsigned *b = LOOP_MALLOC(unsigned , n );
  CUDACALL( cudaMemcpy(b , a , n*sizeof(unsigned) , cudaMemcpyDeviceToDevice));

  return b;
}

unsigned uints_max( unsigned  const* a, unsigned n)
{
  unsigned max = 0;
  thrust::device_ptr< unsigned int const> p (a);
  max = thrust::reduce( p, p +n , INT_MIN , thrust::maximum<unsigned>());
  return max;
}



unsigned* uints_exscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, (n + 1));
  unsigned sum = 0;
  o[0] = 0;
  for (unsigned i = 0; i < n; ++i) {
    sum += a[i];
    o[i + 1] = sum;
  }
  return o;
}
// I hate const
unsigned* uints_unscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  thrust::device_ptr<unsigned> p2 (o);
  thrust::device_ptr<unsigned const> p1 (a);

  thrust::minus<unsigned int> op;

  thrust::transform( p1, p1+n ,p1+1,  p2 , op );


  return o;
}


unsigned* uints_negate(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  thrust::device_ptr<unsigned> p2 (o);
  thrust::device_ptr<unsigned const> p1 (a);

  thrust::transform( p1, p1 +n , p2 , thrust::negate<unsigned>());
  return o;
}


unsigned* uints_negate_offsets(unsigned const* a, unsigned n)
{
  unsigned* unscanned = uints_unscan(a, n);
  unsigned* negated = uints_negate(unscanned, n);
  loop_free(unscanned);
  unsigned* out = uints_exscan(negated, n);
  loop_free(negated);
  return out;
}

void uints_fill(unsigned* a, unsigned n, unsigned v)
{
	thrust::device_ptr< unsigned int> p (a);
	thrust::fill( p , p+n , v);
}

unsigned uints_sum(unsigned const* a, unsigned n)
{
	  unsigned sum = 0;
	  thrust::device_ptr< unsigned int const> p (a);

	  sum = thrust::reduce( p, p +n , (unsigned)0 , thrust::plus<unsigned>());
	  return sum;
}

static int uints_less(void const* a, void const* b)
{
  unsigned const* pa = (unsigned const*) a;
  unsigned const* pb = (unsigned const*) b;
  if (*pa < *pb)
    return -1;
  if (*pa > *pb)
    return 1;
  return 0;
}

unsigned* uints_sort(unsigned const* a, unsigned n)
{
  unsigned* out = uints_copy(a, n);
  thrust::device_ptr< unsigned int> p (out);
  thrust::sort( p , n+p);
  return out;
}

void uints_unique(unsigned const* a, unsigned n,
    unsigned* nunique, unsigned** unique)
{
  unsigned* sorted = uints_sort(a, n);
  unsigned* jump = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    jump[i] = ((i == 0) || (sorted[i - 1] != sorted[i]));
  unsigned* scan = uints_exscan(jump, n);
  *nunique = scan[n];
  *unique = LOOP_MALLOC(unsigned, *nunique);
  for (unsigned i = 0; i < n; ++i)
    if (jump[i])
      (*unique)[scan[i]] = sorted[i];
  loop_free(sorted);
  loop_free(jump);
  loop_free(scan);
}



unsigned long* ulongs_copy(unsigned long const * a , unsigned n)
{
  unsigned long *b = LOOP_MALLOC(unsigned long , n );
  CUDACALL( cudaMemcpy(b, a, n*sizeof(unsigned long), cudaMemcpyDeviceToDevice));
  return b;
}




unsigned long ulongs_max(unsigned long const* a, unsigned n)
{
	  unsigned long max = 0;
	  thrust::device_ptr< unsigned long const> p (a);
	  max = thrust::reduce( p, p +n , LONG_MIN , thrust::maximum<unsigned long>());
	  return max;
}


unsigned char* uchars_copy(unsigned char const* a, unsigned n)
{
	unsigned char *b = LOOP_MALLOC( unsigned char, n);
	CUDACALL( cudaMemcpy(b,a,n*sizeof(unsigned char) , cudaMemcpyDeviceToDevice));
	return b;
}

#else

void uints_zero(unsigned* a, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    a[i] = 0;
}

unsigned* uints_copy(unsigned const* a, unsigned n)
{
  unsigned* b = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}

unsigned* uints_copy(unsigned * a, unsigned n)
{
  unsigned* b = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}

unsigned uints_max(unsigned const* a, unsigned n)
{
  unsigned max = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

unsigned* uints_exscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, (n + 1));
  unsigned sum = 0;
  o[0] = 0;
  for (unsigned i = 0; i < n; ++i) {
    sum += a[i];
    o[i + 1] = sum;
  }
  return o;
}

unsigned* uints_unscan(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    o[i] = a[i + 1] - a[i];
  return o;
}

unsigned* uints_negate(unsigned const* a, unsigned n)
{
  unsigned* o = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    o[i] = !a[i];
  return o;
}

unsigned* uints_negate_offsets(unsigned const* a, unsigned n)
{
  unsigned* unscanned = uints_unscan(a, n);
  unsigned* negated = uints_negate(unscanned, n);
  loop_free(unscanned);
  unsigned* out = uints_exscan(negated, n);
  loop_free(negated);
  return out;
}

void uints_fill(unsigned* a, unsigned n, unsigned v)
{
  for (unsigned i = 0; i < n; ++i)
    a[i] = v;
}

unsigned uints_sum(unsigned const* a, unsigned n)
{
  unsigned sum = 0;
  for (unsigned i = 0; i < n; ++i)
    sum += a[i];
  return sum;
}

static int uints_less(void const* a, void const* b)
{
  unsigned const* pa = (unsigned const*) a;
  unsigned const* pb = (unsigned const*) b;
  if (*pa < *pb)
    return -1;
  if (*pa > *pb)
    return 1;
  return 0;
}

unsigned* uints_sort(unsigned const* a, unsigned n)
{
  unsigned* out = uints_copy(a, n);
  qsort(out, n, sizeof(unsigned), uints_less);
  return out;
}

void uints_unique(unsigned const* a, unsigned n,
    unsigned* nunique, unsigned** unique)
{
  unsigned* sorted = uints_sort(a, n);
  unsigned* jump = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    jump[i] = ((i == 0) || (sorted[i - 1] != sorted[i]));
  unsigned* scan = uints_exscan(jump, n);
  *nunique = scan[n];
  *unique = LOOP_MALLOC(unsigned, *nunique);
  for (unsigned i = 0; i < n; ++i)
    if (jump[i])
      (*unique)[scan[i]] = sorted[i];
  loop_free(sorted);
  loop_free(jump);
  loop_free(scan);
}

unsigned long* ulongs_copy(unsigned long const* a, unsigned n)
{
  unsigned long* b = LOOP_MALLOC(unsigned long, n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}

unsigned long ulongs_max(unsigned long const* a, unsigned n)
{
  unsigned long max = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] > max)
      max = a[i];
  return max;
}

unsigned char* uchars_copy(unsigned char const* a, unsigned n)
{
  unsigned char* b = LOOP_MALLOC(unsigned char, n);
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
  return b;
}
#endif
