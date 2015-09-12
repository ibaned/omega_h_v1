#include "recover_by_volume.h"
#include <string.h>  // for strcat, strcpy
#include <stdlib.h>  // for malloc, free
#include <string.h>  // for strlen
#include "field.h"   // for const_field
#include "mesh.h"    // for mesh_ask_up, mesh_dim, mesh_find_elem_field
#include "size.h"    // for mesh_element_sizes

double* recover_by_volume(
    unsigned nverts,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    double const* size_of_elems,
    unsigned ncomps,
    double const* comps_of_elems)
{
  double* comps_of_verts = malloc(sizeof(double) * ncomps * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
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
  return comps_of_verts;
}

struct const_field* mesh_recover_by_volume(
    struct mesh* m, char const* name)
{
  mesh_element_sizes(m);
  double const* elem_sizes = mesh_find_elem_field(m, "elem_size")->data;
  struct const_field* f = mesh_find_elem_field(m, name);
  double* data = recover_by_volume(mesh_count(m, 0),
      mesh_ask_up(m, 0, mesh_dim(m))->offsets,
      mesh_ask_up(m, 0, mesh_dim(m))->adj,
      elem_sizes,
      f->ncomps, f->data);
  mesh_free_elem_field(m, "elem_size");
  static char const* prefix = "rcov_";
  char* rcov_name = malloc(strlen(name) + strlen(prefix) + 1);
  strcpy(rcov_name, prefix);
  strcat(rcov_name, name);
  struct const_field* out = mesh_add_nodal_field(m, rcov_name, f->ncomps, data);
  free(rcov_name);
  return out;
}
