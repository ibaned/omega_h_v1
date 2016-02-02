#include "gmsh_io.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "files.h"
#include "find_by_verts.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

static unsigned get_gmsh_type_dim(unsigned type)
{
  switch (type) {
    case 15: return 0;
    case 1: return 1;
    case 2: return 2;
    case 4: return 3;
    default: return INVALID;
  }
}

static unsigned find_nodes_start(FILE* f)
{
  static char const* const plain_prefix = "$Nodes";
  static char const* const param_prefix = "$ParametricNodes";
  unsigned long plain_len = strlen(plain_prefix);
  unsigned long param_len = strlen(param_prefix);
  line_t line;
  while (fgets(line, (int) sizeof(line), f)) {
    if (!strncmp(line, plain_prefix, plain_len))
      return 0;
    if (!strncmp(line, param_prefix, param_len))
      return 1;
  }
  LOOP_NORETURN(0);
}

static void read_nodes(FILE* f, unsigned* p_nnodes,
    double** p_node_coords, double** p_node_params)
{
  unsigned has_param = find_nodes_start(f);
  unsigned nnodes;
  safe_scanf(f, 1, "%u", &nnodes);
  double* node_coords = LOOP_HOST_MALLOC(double, nnodes * 3);
  double* node_params = 0;
  if (has_param)
    node_params = LOOP_HOST_MALLOC(double, nnodes * 2);
  for (unsigned i = 0; i < nnodes; ++i) {
    unsigned node_id;
    safe_scanf(f, 4, "%u %lf %lf %lf",
        &node_id,
        node_coords + i * 3 + 0,
        node_coords + i * 3 + 1,
        node_coords + i * 3 + 2);
    assert(node_id == i + 1);
    unsigned geom_dim;
    if (has_param) {
      safe_scanf(f, 1, "%u %*u", &geom_dim);
      switch (geom_dim) {
        case 0:
        case 3: node_params[i * 2 + 0] = 0.0; 
                node_params[i * 2 + 1] = 0.0; 
                break;
        case 1: safe_scanf(f, 1, "%lf", node_params + i * 2);
                node_params[i * 2 + 1] = 0.0;
                break;
        case 2: safe_scanf(f, 2, "%lf %lf",
                    node_params + i * 2 + 0,
                    node_params + i * 2 + 1);
                break;
        default: assert(0);
      }
    }
  }
  *p_nnodes = nnodes;
  *p_node_coords = node_coords;
  *p_node_params = node_params;
}

struct mesh* read_msh(char const* filename)
{
  FILE* f = safe_fopen(filename, "r");
  assert(f != NULL);
  unsigned nnodes;
  double* node_coords;
  double* node_params;
  read_nodes(f, &nnodes, &node_coords, &node_params);
  line_t line;
  seek_prefix(f, line, sizeof(line), "$Elements");
  /* we call these "eq"s, for equal-order classified mesh entities */
  unsigned neqs;
  safe_scanf(f, 1, "%u", &neqs);
  /* conservative allocation, each entity will get 4 vertex
     slots and may leave the later ones unused */
  unsigned* verts_of_eqs = LOOP_HOST_MALLOC(unsigned, neqs * 4);
  unsigned* class_id_of_eqs = LOOP_HOST_MALLOC(unsigned, neqs);
  unsigned* dim_of_eqs = LOOP_HOST_MALLOC(unsigned, neqs);
  unsigned last_dim = 0;
  for (unsigned i = 0; i < neqs; ++i) {
    unsigned type, ntags;
    safe_scanf(f, 2, "%*u %u %u", &type, &ntags);
    dim_of_eqs[i] = get_gmsh_type_dim(type);
    assert(dim_of_eqs[i] >= last_dim);
    assert(dim_of_eqs[i] <= last_dim + 1);
    assert(ntags >= 2);
    unsigned physical;
    /* grab the physical and raw geometric classifications */
    safe_scanf(f, 2, "%u %u", &physical, class_id_of_eqs + i);
    /* for now ignore physical, but will be very useful in the future */
    (void) physical;
    for (unsigned j = 2; j < ntags; ++j)
      safe_scanf(f, 0, "%*u"); /* discard other "tags" */
    unsigned nent_verts = dim_of_eqs[i] + 1;
    for (unsigned j = 0; j < nent_verts; ++j) {
      safe_scanf(f, 1, "%u", verts_of_eqs + i * 4 + j);
      verts_of_eqs[i * 4 + j] -= 1;
    }
    last_dim = dim_of_eqs[i];
  }
  fclose(f);
  /* alright, we can tell the highest-dimensional entity
     and store the vertices at least */
  unsigned dim = uints_max(dim_of_eqs, neqs);
  struct mesh* m = new_mesh(dim, MESH_REDUCED, 0);
  mesh_set_ents(m, 0, nnodes, 0);
  mesh_add_tag(m, 0, TAG_F64, "coordinates", 3, node_coords);
  if (node_params)
    mesh_add_tag(m, 0, TAG_F64, "parametric", 2, node_params);
  /* we can also form the elements and derive all intermediate
     entities based on the elements */
  unsigned nelems = 0;
  for (unsigned i = 0; i < neqs; ++i)
    if (dim_of_eqs[i] == dim)
      ++nelems;
  unsigned verts_per_elem = dim + 1;
  unsigned* verts_of_elems = LOOP_HOST_MALLOC(unsigned, nelems * verts_per_elem);
  unsigned ei = 0;
  for (unsigned i = 0; i < neqs; ++i) {
    if (dim_of_eqs[i] != dim)
      continue;
    for (unsigned j = 0; j < verts_per_elem; ++j)
      verts_of_elems[ei * verts_per_elem + j] =
          verts_of_eqs[i * 4 + j];
    ++ei;
  }
  mesh_set_ents(m, dim, nelems, verts_of_elems);
  /* now for the tricky bit. for each equal-order entity from the file,
     find its mesh-structure counterpart, and classify the closure of
     that mesh entity the same as the equal-order classification.
     to get the right results, we have to do this in order of decreasing
     dimension, which as we've asserted while reading the file, is
     decreasing order of appearance in the file */
  unsigned* class_dims[4];
  unsigned* class_ids[4];
  for (unsigned i = 0; i <= dim; ++i) {
    unsigned n = mesh_count(m, i);
    class_dims[i] = LOOP_MALLOC(unsigned, n);
    class_ids[i] = LOOP_MALLOC(unsigned, n);
  }
  for (unsigned ii = 0; ii < neqs; ++ii) {
    unsigned i = neqs - ii - 1;
    unsigned eq_dim = dim_of_eqs[i];
    unsigned const* eq_verts = verts_of_eqs + i * 4;
    /* general reverse lookup via vertices. */
    unsigned ent = find_by_verts(eq_dim + 1, eq_verts,
        mesh_ask_down(m, eq_dim, 0),
        mesh_ask_up(m, 0, eq_dim)->adj,
        mesh_ask_up(m, 0, eq_dim)->offsets);
    assert(ent != INVALID);
    class_dims[eq_dim][ent] = eq_dim;
    class_ids[eq_dim][ent] = class_id_of_eqs[i];
    for (unsigned dd = 0; dd < eq_dim; ++dd) {
      unsigned const* des_of_ents = mesh_ask_down(m, eq_dim, dd);
      unsigned des_per_ent = the_down_degrees[eq_dim][dd];
      for (unsigned j = 0; j < des_per_ent; ++j) {
        unsigned de = des_of_ents[ent * des_per_ent + j];
        class_dims[dd][de] = eq_dim;
        class_ids[dd][de] = class_id_of_eqs[i];
      }
    }
  }
  loop_free(verts_of_eqs);
  loop_free(class_id_of_eqs);
  loop_free(dim_of_eqs);
  for (unsigned i = 0; i <= dim; ++i) {
    mesh_add_tag(m, i, TAG_U32, "class_dim", 1, class_dims[i]);
    mesh_add_tag(m, i, TAG_U32, "class_id", 1, class_ids[i]);
  }
  mesh_set_rep(m, MESH_FULL);
  return m;
}
