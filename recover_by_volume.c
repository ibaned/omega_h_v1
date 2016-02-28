#include "recover_by_volume.h"

#include "loop.h"
#include "mesh.h"
#include "size.h"
#include "tag.h"

LOOP_KERNEL(execute,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    double const* size_of_elems,
    unsigned ncomps,
    double const* comps_of_elems,
    double* comps_of_verts)
  
  unsigned first_use = elems_of_verts_offsets[i];
  unsigned end_use = elems_of_verts_offsets[i + 1];
  double* comps_of_vert = comps_of_verts + i * ncomps;
  for (unsigned j = 0; j < ncomps; ++j)
    comps_of_vert[j] = 0;
  double size_sum = 0;
  for (unsigned j = first_use; j < end_use; ++j) {
    unsigned elem = elems_of_verts[j];
    double elem_size = size_of_elems[elem];
    size_sum += elem_size;
    double const* comps_of_elem = comps_of_elems + elem * ncomps;
    for (unsigned k = 0; k < ncomps; ++k)
      comps_of_vert[k] += elem_size * comps_of_elem[k];
  }
  for (unsigned j = 0; j < ncomps; ++j)
    comps_of_vert[j] /= size_sum;
}

double* recover_by_volume(
    unsigned nverts,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    double const* size_of_elems,
    unsigned ncomps,
    double const* comps_of_elems)
{
  double* comps_of_verts = LOOP_MALLOC(double, ncomps * nverts);
  LOOP_EXEC( execute, nverts,
    elems_of_verts_offsets,
    elems_of_verts,
    size_of_elems,
    ncomps,
    comps_of_elems,
    comps_of_verts);
  return comps_of_verts;
}

struct const_tag* mesh_recover_by_volume(
    struct mesh* m, char const* name)
{
  double* elem_sizes = mesh_element_sizes(m);
  struct const_tag* t = mesh_find_tag(m, mesh_dim(m), name);
  double* data = recover_by_volume(mesh_count(m, 0),
      mesh_ask_up(m, 0, mesh_dim(m))->offsets,
      mesh_ask_up(m, 0, mesh_dim(m))->adj,
      elem_sizes,
      t->ncomps, t->d.f64);
  loop_free(elem_sizes);
  struct const_tag* out = mesh_add_tag(m, 0, TAG_F64, name, t->ncomps, data);
  return out;
}
