#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "classify_box.h"
#include "eval_field.h"
#include "mesh.h"
#include "refine_by_size.h"
#include "vtk.h"

static double size = 0.75;

static void sizefun(double const* x, double* o)
{
  (void) x;
  *o = size;
}

int main(int argc, char** argv)
{
  unsigned dim = 3;
  char const* file = "out.vtu";
  for (int i = 0; i < argc; ++i) {
    if (!strcmp(argv[i], "--dim")) {
      ++i;
      if (i == argc) {
        printf("--dim needs an argument\n");
        return -1;
      }
      dim = (unsigned) atoi(argv[i]);
      assert(dim <= 3);
    } else if (!strcmp(argv[i], "--size")) {
      ++i;
      if (i == argc) {
        printf("--size needs an argument\n");
        return -1;
      }
      size = atof(argv[i]);
    } else if (!strcmp(argv[i], "--file")) {
      ++i;
      if (i == argc) {
        printf("--file needs an argument\n");
        return -1;
      }
      file = argv[i];
    }
  }
  struct mesh* m = new_box_mesh(dim);
  mesh_eval_field(m, 0, "adapt_size", 1, sizefun);
  while (refine_by_size(&m, 0.0));
  mesh_free_tag(m, 0, "adapt_size");
  mesh_classify_box(m);
  write_vtu(m, file);
  free_mesh(m);
}
