#include "element_field.h"

#include "loop.h"
#include "mesh.h"
#include "tables.h"
#include "tag.h"

LOOP_KERNEL(interp_to_elem,
    unsigned const* verts_of_elems,
    unsigned verts_per_elem,
    double const* in,
    unsigned ncomps,
    double* out)
  unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
  average_element_field(verts_per_elem, verts_of_elem, in, ncomps,
      out + i * ncomps);
}

double* interp_to_elems(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned ncomps,
    double const* in)
{
  double* out = LOOP_MALLOC(double, nelems * ncomps);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  LOOP_EXEC(interp_to_elem, nelems, verts_of_elems, verts_per_elem,
      in, ncomps, out);
  return out;
}

void mesh_interp_to_elems(struct mesh* m, char const* name)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  struct const_tag* t = mesh_find_tag(m, 0, name);
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  double* data = interp_to_elems(elem_dim, nelems, verts_of_elems,
      t->ncomps, t->d.f64);
  mesh_add_tag(m, elem_dim, TAG_F64, t->name, t->ncomps, data);
}
