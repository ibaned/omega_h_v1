#include "vtk.h"
#include "mesh.h"
#include "mark.h"
#include "subset.h"
#include "loop.h"
#include "ints.h"
#include "field.h"
#include "recover_by_volume.h"
#include <assert.h>
#include <stdio.h>

static void list_fields(struct mesh* m)
{
  for (unsigned i = 0; i < mesh_count_nodal_fields(m); ++i)
    printf("nodal field \"%s\"\n",
        mesh_get_nodal_field(m, i)->name);
  for (unsigned i = 0; i < mesh_count_elem_fields(m); ++i)
    printf("element field \"%s\"\n",
        mesh_get_elem_field(m, i)->name);
}

int main(int argc, char** argv)
{
  assert(argc == 3);
  struct mesh* m = read_vtk(argv[1]);
  for (unsigned i = 0; i < mesh_count_elem_fields(m); ++i) {
    struct const_field* f = mesh_get_elem_field(m, i);
    mesh_recover_by_volume(m, f->name);
  }
  list_fields(m);
  unsigned d = mesh_dim(m);
  unsigned ns = mesh_count(m, d - 1);
  unsigned* pb = mesh_mark_part_boundary(m);
  unsigned* off = ints_exscan(pb, ns);
  loop_free(pb);
  struct mesh* sm = subset_mesh(m, d - 1, off);
  list_fields(sm);
  free_mesh(m);
  loop_free(off);
  write_vtk(sm, argv[2]);
  free_mesh(sm);
}
