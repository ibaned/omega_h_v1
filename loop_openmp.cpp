#include <omp.h>

unsigned loop_size(void)
{
  return (unsigned) omp_get_max_threads();
}
