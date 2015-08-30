#include "subset.h"
#include "mesh.h"
#include "tables.h"
#include "ints.h"
#include <stdlib.h>
#include <string.h>

static void* general_subset(
    unsigned typesize,
    unsigned n,
    unsigned width,
    void const* a,
    unsigned const* offsets)
{
  unsigned nsub = offsets[n];
  void* out = malloc(typesize * nsub * width);
  unsigned stride = typesize * width;
  char const* ip = a;
  char* op = out;
  for (unsigned i = 0; i < n; ++i, ip += stride) {
    if (offsets[i] == offsets[i + 1])
      continue;
    memcpy(op, ip, stride);
    op += stride;
  }
  return out;
}

unsigned* ints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* offsets)
{
  return general_subset(sizeof(unsigned), n, width, a, offsets);
}

double* doubles_subset(
    unsigned n,
    unsigned width,
    double const* a,
    unsigned const* offsets)
{
  return general_subset(sizeof(double), n, width, a, offsets);
}

static unsigned* mark_elem_verts(
    unsigned nverts,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    unsigned const* offsets)
{
  unsigned* marked = malloc(sizeof(unsigned) * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_use = elems_of_verts_offsets[i];
    unsigned end_use = elems_of_verts_offsets[i + 1];
    marked[i] = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned elem = elems_of_verts[j];
      if (offsets[elem] != offsets[elem + 1]) {
        marked[i] = 1;
        break;
      }
    }
  }
  unsigned* out = ints_exscan(marked, nverts);
  free(marked);
  return out;
}

struct mesh* subset_mesh(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* offsets)
{
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned nelems_out = offsets[nelems];
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* verts_of_elems_out = ints_subset(nelems, verts_per_elem,
      verts_of_elems, offsets);
  unsigned nverts = mesh_count(m, 0);
  unsigned* vert_offsets = mark_elem_verts(nverts,
      mesh_ask_up(m, 0, elem_dim)->offsets,
      mesh_ask_up(m, 0, elem_dim)->adj,
      offsets);
  for (unsigned i = 0; i < nelems_out * verts_per_elem; ++i)
    verts_of_elems_out[i] = vert_offsets[verts_of_elems_out[i]];
  unsigned nverts_out = vert_offsets[nverts];
  struct mesh* out = new_mesh(elem_dim);
  mesh_set_ents(out, 0, nverts_out, 0);
  mesh_set_ents(out, elem_dim, nelems_out, verts_of_elems_out);
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  double* coords_out = doubles_subset(nverts, 3, coords, vert_offsets);
  mesh_add_nodal_field(out, "coordinates", 3, coords_out);
  return out;
}
