#include "gmsh_io.h"

#include <assert.h>
#include <stdio.h>

#include "files.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

static unsigned get_dim_type(unsigned dim)
{
  switch (dim) {
    case 0: return 15;
    case 1: return 1;
    case 2: return 2;
    case 3: return 4;
    default: return 0;
  }
}

struct mesh* read_msh(char const* filename)
{
  FILE* f = fopen(filename, "r");
  assert(f != NULL);
  line_t line;
  seek_prefix(f, line, sizeof(line), "$Nodes");
  unsigned nnodes;
  safe_scanf(f, 1, "%u", &nnodes);
  double* node_coords = LOOP_HOST_MALLOC(double, nnodes * 3);
  for (unsigned i = 0; i < nnodes; ++i) {
    unsigned node_id;
    safe_scanf(f, 4, "%u %lf %lf %lf",
        &node_id,
        node_coords + i * 3 + 0,
        node_coords + i * 3 + 1,
        node_coords + i * 3 + 2);
    assert(node_id == i + 1);
  }
  unsigned dim = 0;
  for (unsigned j = 0; j < 3; ++j) {
    for (unsigned i = 0; i < nnodes; ++i) {
      if (node_coords[i * 3 + j] != 0.0) {
        ++dim;
        break;
      }
    }
  }
  assert(dim > 0);
  assert(dim <= 3);
  unsigned* class_dim = LOOP_HOST_MALLOC(unsigned, nnodes);
  unsigned* class_id = LOOP_HOST_MALLOC(unsigned, nnodes);
  unsigned* class_phys = LOOP_HOST_MALLOC(unsigned, nnodes);
  for (unsigned i = 0; i < nnodes; ++i)
    class_dim[i] = INVALID;
  seek_prefix(f, line, sizeof(line), "$Elements");
  unsigned nents;
  safe_scanf(f, 1, "%u", &nents);
  unsigned nverts_per_elem = dim + 1;
  /* conservative allocation, nents >= nelems */
  unsigned* verts_of_elems = LOOP_HOST_MALLOC(unsigned, nents * nverts_per_elem);
  unsigned ent_dim = 0;
  unsigned nelems = 0;
  for (unsigned i = 0; i < nents; ++i) {
    unsigned type, ntags;
    safe_scanf(f, 2, "%*u %u %u", &type, &ntags);
    if (type != get_dim_type(ent_dim)) {
      ++ent_dim;
      assert(type == get_dim_type(ent_dim));
    }
    assert(ntags >= 2);
    unsigned physical, cad;
    safe_scanf(f, 2, "%u %u", &physical, &cad);
    for (unsigned j = 2; j < ntags; ++j)
      safe_scanf(f, 0, "%*u");
    unsigned nent_verts = ent_dim + 1;
    unsigned ent_verts[4];
    for (unsigned j = 0; j < nent_verts; ++j) {
      safe_scanf(f, 1, "%u", ent_verts + j);
      ent_verts[j] -= 1;
    }
    for (unsigned j = 0; j < nent_verts; ++j)
      if (class_dim[ent_verts[j]] == INVALID) {
        class_dim[ent_verts[j]] = ent_dim;
        class_id[ent_verts[j]] = cad;
        class_phys[ent_verts[j]] = physical;
      }
    if (ent_dim == dim) {
      for (unsigned j = 0; j < nent_verts; ++j)
        verts_of_elems[nelems * nverts_per_elem + j] =
            ent_verts[j];
      ++nelems;
    }
  }
  fclose(f);
  struct mesh* m = new_mesh(dim);
  mesh_set_ents(m, 0, nnodes, 0);
  mesh_add_tag(m, 0, TAG_F64, "coordinates", 3, node_coords);
  mesh_add_tag(m, 0, TAG_U32, "class_dim", 1, class_dim);
  mesh_add_tag(m, 0, TAG_U32, "class_id", 1, class_id);
  mesh_add_tag(m, 0, TAG_U32, "class_phys", 1, class_phys);
  mesh_set_ents(m, dim, nelems, verts_of_elems);
  return m;
}
