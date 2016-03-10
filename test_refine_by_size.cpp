#include <math.h>
#include <stdio.h>

#include "algebra.h"
#include "arrays.h"
#include "comm.h"
#include "derive_model.h"
#include "doubles.h"
#include "eval_field.h"
#include "loop.h"
#include "mesh.h"
#include "quality.h"
#include "refine.h"
#include "tag.h"
#include "vtk_io.h"

static void size_fun(double const* x, double* s)
{
  double coarse = 0.5;
  double fine = 0.1;
  double radius = vector_norm(x, 3);
  double d = fabs(radius - 0.5);
  s[0] = coarse * d + fine * (1 - d);
}

int main()
{
  comm_init();
  struct mesh* m = new_box_mesh(2);
  mesh_derive_model(m, PI / 4);
  mesh_set_rep(m, MESH_FULL);
  char fname[64];
  mesh_eval_field(m, 0, "adapt_size", 1, size_fun);
  { //set mass field to test conservative transfer
    unsigned nelems = mesh_count(m, mesh_dim(m));
    mesh_add_tag(m, mesh_dim(m), TAG_F64, "mass", 1,
        doubles_filled(nelems, 1.0 / nelems));
  }
  write_mesh_vtk(m, "out_0.vtu");
  for (unsigned it = 1; 1; ++it) {
    if (!refine_by_size(m, 0))
      break;
    printf("%u elements, %u vertices\n", mesh_count(m, mesh_dim(m)), mesh_count(m, 0));
    double* quals = element_qualities(mesh_dim(m), mesh_count(m, mesh_dim(m)),
        mesh_ask_down(m, mesh_dim(m), 0),
        mesh_find_tag(m, 0, "coordinates")->d.f64);
    double minqual = doubles_min(quals, mesh_count(m, mesh_dim(m)));
    loop_free(quals);
    printf("min quality %f\n", minqual);
    sprintf(fname, "out_%u.vtu", it);
    write_mesh_vtk(m, fname);
    mesh_free_tag(m, 0, "adapt_size");
    mesh_eval_field(m, 0, "adapt_size", 1, size_fun);
  }
  free_mesh(m);
  comm_fini();
}
