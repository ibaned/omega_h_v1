#include "include/omega_h.h"
#include "comm.hpp"

#include <cstdio>
#include <cassert>

using namespace omega_h;

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  unsigned world_size = comm_size();
  unsigned subcomm_size = 0;
  for (subcomm_size = 1; subcomm_size <= world_size; subcomm_size *= 2) {
    unsigned group = comm_rank() / subcomm_size;
    struct comm* subcomm = comm_split(comm_using(), group, comm_rank() % subcomm_size);
    comm_use(subcomm);
    if (!group)
      printf("am now subcomm rank %u of %u\n", comm_rank(), comm_size());
    comm_use(comm_world());
    comm_free(subcomm);
  }
  osh_fini();
}
