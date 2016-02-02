#include <math.h>
#include <stdio.h>
//#include <mpi.h>

#include "comm.h"
#include "eval_field.h"
#include "mesh.h"
#include "refine.h"
#include "vtk.h"

static void sinusoid(double const* x, double* size)
{
  double s = cos(x[0] * 8.0 * M_PI) / 4.0 + 1.0 / 2.0;
  double d = fabs(x[1] - s);
  *size = sqrt(2 * (1e-6 + d * 1e-5));
}

int main()
{
  comm_init();
  struct mesh* m = new_box_mesh(2);
  unsigned did_refine = 0;
  //double t0 = MPI_Wtime();
  do {
    mesh_eval_field(m, 0, "adapt_size", 1, sinusoid);
    did_refine = refine_by_size(&m, 0.0);
    mesh_free_tag(m, 0, "adapt_size");
  } while (did_refine);
  //double t1 = MPI_Wtime();
  //printf("refinement time %.3e seconds, %u final triangles\n",
  //    t1 - t0, mesh_count(m, 2));
  write_mesh_vtk(m, "flounder.vtu");
  free_mesh(m);
  comm_fini();
}
