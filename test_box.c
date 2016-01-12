#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "comm.h"
#include "eval_field.h"
#include "mesh.h"
#include "refine_by_size.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  comm_init();
  unsigned dim = 3;
  unsigned nrefs = 0;
  char const* file = "out.vtu";
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--dim")) {
      ++i;
      if (i == argc) {
        printf("--dim needs an argument\n");
        return -1;
      }
      dim = (unsigned) atoi(argv[i]);
      assert(dim <= 3);
    } else if (!strcmp(argv[i], "--refinements")) {
      ++i;
      if (i == argc) {
        printf("--refinements needs an argument\n");
        return -1;
      }
      nrefs = (unsigned) atoi(argv[i]);
      assert(nrefs < 30);
    } else if (!strcmp(argv[i], "--file")) {
      ++i;
      if (i == argc) {
        printf("--file needs an argument\n");
        return -1;
      }
      file = argv[i];
    } else {
      printf("unknown argument %s\n", argv[i]);
      return -1;
    }
  }
  struct mesh* m = new_box_mesh(dim);
  for (unsigned i = 0; i < nrefs; ++i)
    uniformly_refine(&m);
  write_mesh_vtk(m, file);
  free_mesh(m);
  comm_fini();
}
