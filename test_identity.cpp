#include <assert.h>
#include <stdlib.h>

#include "include/omega_h.hpp"

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  assert(argc == 3);
  osh_t m = osh_read_vtk(argv[1]);
  osh_identity_size(m, "adapt_size");
  osh_adapt(m, 1./3., 0.3, 4, 10);
  osh_write_vtk(m, argv[2]);
  osh_free(m);
  osh_fini();
}

