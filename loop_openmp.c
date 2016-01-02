#include <omp.h>

unsigned loop_size(void)
{
  return omp_get_max_threads();
}
