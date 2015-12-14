#include "include/omega_h.h"

#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "parallel_mesh.h"
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

osh_t osh_build(unsigned dim, unsigned nelems, unsigned nverts,
    unsigned* conn)
{
  struct mesh* m = new_mesh(dim);
  mesh_set_ents(m, 0, nverts, 0);
  mesh_set_ents(m, 0, nelems, conn);
  return (osh_t) m;
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

unsigned const* osh_own_rank(osh_t m, unsigned dim)
{
  return mesh_ask_own_ranks((struct mesh*)m, dim);
}

unsigned const* osh_own_id(osh_t m, unsigned dim)
{
  return mesh_ask_own_ids((struct mesh*)m, dim);
}

void osh_set_field(osh_t m, char const* name, unsigned ncomps, double* data)
{
  struct const_tag* t = mesh_find_tag((struct mesh*)m, 0, name);
  if (t)
    mesh_free_tag((struct mesh*)m, 0, name);
  mesh_add_tag((struct mesh*)m, 0, TAG_F64, name, ncomps, data);
}

void osh_new_field(osh_t m, char const* name, unsigned ncomps)
{
  if (mesh_find_tag((struct mesh*)m, 0, name))
    return;
  double* data = LOOP_MALLOC(double, ncomps * mesh_count((struct mesh*)m, 0));
  mesh_add_tag((struct mesh*)m, 0, TAG_F64, name, ncomps, data);
}

double* osh_get_field(osh_t m, char const* name)
{
  return mesh_find_tag((struct mesh*)m, 0, name)->d.f64;
}

void osh_new_label(osh_t m, char const* name, unsigned ncomps)
{
  if (mesh_find_tag((struct mesh*)m, 0, name))
    return;
  unsigned* data = LOOP_MALLOC(unsigned, ncomps * mesh_count((struct mesh*)m, 0));
  mesh_add_tag((struct mesh*)m, 0, TAG_U32, name, ncomps, data);
}

unsigned* osh_get_label(osh_t m, char const* name)
{
  return mesh_find_tag((struct mesh*)m, 0, name)->d.u32;
}

void osh_set_global(osh_t m, unsigned long* data)
{
  mesh_add_tag((struct mesh*)m, 0, TAG_U64, "global_number", 1, data);
}

void osh_mark_verts(osh_t m, unsigned class_dim, unsigned class_id,
    unsigned* marked)
{
  unsigned* to_mark = mesh_mark_class((struct mesh*)m, 0,
      class_dim, class_id);
  unsigned nverts = osh_nverts(m);
  for (unsigned i = 0; i < nverts; ++i)
    if (to_mark[i])
      marked[i] = 1;
  loop_free(to_mark);
}

void osh_add_label(osh_t m, char const* name, unsigned* data)
{
  mesh_add_tag((struct mesh*)m, 0, TAG_U32, name, 1, data);
}

void osh_free_label(osh_t m, char const* name)
{
  mesh_free_tag((struct mesh*)m, 0, name);
}
