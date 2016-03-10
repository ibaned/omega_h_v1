#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "algebra.hpp"
#include "derive_model.hpp"
#include "eval_field.hpp"
#include "int_casts.hpp"
#include "mesh.hpp"
#include "include/omega_h.hpp"
#include "refine.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  unsigned dim = 3;
  unsigned nrefs = 0;
  unsigned full = 1;
  char const* file = "box.vtu";
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--dim")) {
      ++i;
      if (i == argc) {
        printf("--dim needs an argument\n");
        return -1;
      }
      dim = U(atoi(argv[i]));
      assert(dim <= 3);
    } else if (!strcmp(argv[i], "--refinements")) {
      ++i;
      if (i == argc) {
        printf("--refinements needs an argument\n");
        return -1;
      }
      nrefs = U(atoi(argv[i]));
      assert(nrefs < 30);
    } else if (!strcmp(argv[i], "--file")) {
      ++i;
      if (i == argc) {
        printf("--file needs an argument\n");
        return -1;
      }
      file = argv[i];
    } else if (!strcmp(argv[i], "--reduced")) {
      full = 0;
    } else {
      printf("unknown argument %s\n", argv[i]);
      return -1;
    }
  }
  struct mesh* m = new_box_mesh(dim);
  if (full) {
    mesh_derive_model(m, PI / 4);
    mesh_set_rep(m, MESH_FULL);
  }
  for (unsigned i = 0; i < nrefs; ++i)
    uniformly_refine(m);
  write_mesh_vtk(m, file);
  free_mesh(m);
  osh_fini();
}
