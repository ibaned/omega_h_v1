#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "classify_box.h"
#include "eval_field.h"
#include "mesh.h"
#include "refine_by_size.h"
#include "vtk.h"

int main(int argc, char** argv)
{
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
  mesh_classify_box(m);
  write_vtu(m, file);
  free_mesh(m);
}
