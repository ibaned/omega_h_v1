#include <assert.h>
#include <stdio.h>

#include "comm.h"
#include "loop.h"
#include "parallel_inertial_bisect.h"

int main()
{
  comm_init();
  assert(comm_size() == 2);
  unsigned const n_in[2] = {2,2};
  double const coords_in[2][6] = {
    {0,0,0, 2,0,0},
    {1,0,0, 3,0,0}};
  double const masses_in[2][2] = {
    {1, 1},
    {1, 1}};
  unsigned const orig_ranks_in[2][2] = {
    {0, 0},
    {1, 1}};
  unsigned const orig_ids_in[2][2] = {
    {0, 1},
    {0, 1}};
  unsigned n = n_in[comm_rank()];
  double* coords = loop_to_device(coords_in[comm_rank()],
      sizeof(double) * 3 * n);
  double* masses = loop_to_device(masses_in[comm_rank()],
      sizeof(double) * 3 * n);
  unsigned* orig_ranks = loop_to_device(orig_ranks_in[comm_rank()],
      sizeof(unsigned) * n);
  unsigned* orig_ids = loop_to_device(orig_ids_in[comm_rank()],
      sizeof(unsigned) * n);
  recursive_inertial_bisect(&n, &coords, &masses, &orig_ranks, &orig_ids);
  printf("n %u\n", n);
  for (unsigned i = 0; i < n; ++i)
    printf("%u %u : (%f %f %f) %f %u %u\n",
        comm_rank(), i,
        coords[i * 3 + 0],
        coords[i * 3 + 1],
        coords[i * 3 + 2],
        masses[i],
        orig_ranks[i],
        orig_ids[i]);
  loop_free(coords);
  loop_free(masses);
  loop_free(orig_ranks);
  loop_free(orig_ids);
  comm_fini();
}
