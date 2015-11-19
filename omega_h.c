#include "include/omega_h.h"

#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "tag.h"
#include "vtk.h"

void osh_free(osh_t m)
{
  free_mesh((struct mesh*)m);
}

osh_t osh_read_vtk(char const* filename)
{
  return (osh_t) read_vtu(filename);
}

void osh_write_vtk(osh_t m, char const* filename)
{
  write_vtu((struct mesh*)m, filename);
}

unsigned osh_dim(osh_t m)
{
  return mesh_dim((struct mesh*)m);
}

unsigned osh_nelems(osh_t m)
{
  return mesh_count((struct mesh*)m,
      mesh_dim((struct mesh*)m));
}

unsigned osh_nverts(osh_t m)
{
  return mesh_count((struct mesh*)m, 0);
}

unsigned const* osh_down(osh_t m, unsigned high_dim, unsigned low_dim)
{
  return mesh_ask_down((struct mesh*)m, high_dim, low_dim);
}

unsigned const* osh_up(osh_t m, unsigned low_dim, unsigned high_dim)
{
  return mesh_ask_up((struct mesh*)m, low_dim, high_dim)->adj;
}

unsigned const* osh_up_offs(osh_t m, unsigned low_dim, unsigned high_dim)
{
  return mesh_ask_up((struct mesh*)m, low_dim, high_dim)->offsets;
}

unsigned const* osh_up_dirs(osh_t m, unsigned low_dim, unsigned high_dim)
{
  return mesh_ask_up((struct mesh*)m, low_dim, high_dim)->directions;
}

double const* osh_coords(osh_t m)
{
  return mesh_find_tag((struct mesh*)m, 0, "coordinates")->d.f64;
}

void osh_set_field(osh_t m, char const* name, unsigned ncomps, double* data)
{
  struct const_tag* t = mesh_find_tag((struct mesh*)m, 0, name);
  if (t)
    mesh_free_tag((struct mesh*)m, 0, name);
  mesh_add_tag((struct mesh*)m, 0, TAG_F64, name, ncomps, data);
}

unsigned const* osh_ask_label(osh_t m, char const* name)
{
  return mesh_find_tag((struct mesh*)m, 0, name)->d.u32;
}

void osh_mark_verts(osh_t m, unsigned class_dim, unsigned class_id,
    unsigned* marked)
{
  unsigned* to_mark = mesh_mark_class_closure_verts((struct mesh*)m,
      class_dim, class_id);
  unsigned nverts = osh_nverts(m);
  for (unsigned i = 0; i < nverts; ++i)
    if (to_mark[i])
      marked[i] = 1;
  loop_free(to_mark);
}
