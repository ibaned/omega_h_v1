#include <assert.h>
#include <stdlib.h>

#include "infer_class.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "subset.h"
#include "vtk.h"

int main(int argc, char** argv)
{
  assert(argc == 3);
  const char* infile = argv[1];
  unsigned dim = (unsigned) atoi(argv[2]);
  assert(dim <= 3);
  struct mesh* m = read_vtu(infile);
  mesh_ask_class_dim(m, dim);
  unsigned nents = mesh_count(m, dim);
  unsigned* marked = uints_filled(nents, 1);
  unsigned* offsets = uints_exscan(marked, nents);
  loop_free(marked);
  struct mesh* sm = subset_mesh(m, dim, offsets);
  free_mesh(m);
  loop_free(offsets);
  write_vtu(sm, "inferred.vtu");
  free_mesh(sm);
}
