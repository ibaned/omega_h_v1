#include <math.h>
#include <stdio.h>

#include "algebra.h"
#include "doubles.h"
#include "eval_field.h"
#include "loop.h"
#include "mesh.h"
#include "quality.h"
#include "refine_by_size.h"
#include "tag.h"
#include "vtk.h"

static void size_fun(double const* x, double* s)
{
  double coarse = 0.5;
  double fine = 0.025;
  double radius = vector_norm(x, 3);
  double d = fabs(radius - 0.5);
  s[0] = coarse * d + fine * (1 - d);
}

int main()
{
  struct mesh* m = new_box_mesh(3);
  char fname[64];
  mesh_eval_field(m, 0, "adapt_size", 1, size_fun);
  for (unsigned it = 0; 1; ++it) {
    if (!refine_by_size(&m, 0))
      break;
    printf("%u elements, %u vertices\n", mesh_count(m, mesh_dim(m)), mesh_count(m, 0));
    double* quals = element_qualities(mesh_dim(m), mesh_count(m, mesh_dim(m)),
        mesh_ask_down(m, mesh_dim(m), 0),
        mesh_find_tag(m, 0, "coordinates")->data);
    double minqual = doubles_min(quals, mesh_count(m, mesh_dim(m)));
    LOOP_FREE(quals);
    printf("min quality %f\n", minqual);
    sprintf(fname, "out_%u.vtu", it);
    write_vtu(m, fname);
    mesh_free_tag(m, 0, "adapt_size");
    mesh_eval_field(m, 0, "adapt_size", 1, size_fun);
  }
  free_mesh(m);
}
