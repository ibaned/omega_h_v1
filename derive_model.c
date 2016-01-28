#include "derive_model.h"

#include <assert.h>
#include <math.h>

#include "algebra.h"
#include "arrays.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "tables.h"

/* a mesh will have these dimensions:

   element
   side = element - 1
   hinge = side - 1
*/

LOOP_KERNEL(boundary_edge_normal, 
    unsigned const* boundary_edges,
    unsigned const* verts_of_edges,
    double const* coords,
    double* normals)

  if (!boundary_edges[i])
    return;
  subtract_vectors(coords + verts_of_edges[i * 2 + 1] * 3,
                   coords + verts_of_edges[i * 2 + 0] * 3,
                   normals + i * 3,
                   3);
  normalize_vector(normals + i * 3, normals + i * 3, 3);
}

LOOP_KERNEL(boundary_tri_normal, 
    unsigned const* boundary_tris,
    unsigned const* verts_of_tris,
    double const* coords,
    double* normals)

  if (!boundary_tris[i])
    return;
  double basis[2][3];
  subtract_vectors(coords + verts_of_tris[i * 3 + 1] * 3,
                   coords + verts_of_tris[i * 3 + 0] * 3,
                   basis[0],
                   3);
  subtract_vectors(coords + verts_of_tris[i * 3 + 2] * 3,
                   coords + verts_of_tris[i * 3 + 0] * 3,
                   basis[1],
                   3);
  cross_product(basis[0], basis[1], normals + i * 3);
  normalize_vector(normals + i * 3, normals + i * 3, 3);
}

static double* get_boundary_side_normals(
    unsigned side_dim,
    unsigned nsides,
    unsigned const* boundary_sides,
    unsigned const* verts_of_sides,
    double const* coords)
{
  double* normals = LOOP_MALLOC(double, 3 * nsides);
  switch (side_dim) {
    case 1:
      LOOP_EXEC(boundary_edge_normal, nsides, boundary_sides, verts_of_sides,
          coords, normals);
      break;
    case 2:
      LOOP_EXEC(boundary_tri_normal, nsides, boundary_sides, verts_of_sides,
          coords, normals);
      break;
  }
  return normals;
}

static double* get_hinge_angles(
    unsigned nhinges,
    unsigned const* sides_of_hinges,
    unsigned const* sides_of_hinges_offsets,
    unsigned const* boundary_sides,
    unsigned const* boundary_hinges,
    double const* side_normals)
{
  double* angles = LOOP_MALLOC(double, nhinges);
  for (unsigned i = 0; i < nhinges; ++i) {
    if (!boundary_hinges[i]) {
      angles[i] = 0;
      continue;
    }
    unsigned f = sides_of_hinges_offsets[i];
    unsigned e = sides_of_hinges_offsets[i + 1];
    unsigned nadj = 0;
    unsigned adj[2];
    for (unsigned j = f; j < e; ++j) {
      unsigned side = sides_of_hinges[j];
      if (boundary_sides[side]) {
        assert(nadj < 2);
        adj[nadj] = side;
        ++nadj;
      }
    }
    double d = dot_product(
        side_normals + adj[0] * 3,
        side_normals + adj[1] * 3,
        3);
    angles[i] = acos(d);
  }
  return angles;
}

static unsigned* mark_creases(
    unsigned nhinges,
    double const* hingle_angles,
    double crease_angle)
{
  unsigned* creases = LOOP_MALLOC(unsigned, nhinges);
  for (unsigned i = 0; i < nhinges; ++i)
    creases[i] = hingle_angles[i] > crease_angle;
  return creases;
}

static unsigned* mark_corners(
    unsigned nverts,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_offsets,
    unsigned const* crease_edges)
{
  unsigned* corners = LOOP_MALLOC(unsigned, nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned f = edges_of_verts_offsets[i];
    unsigned e = edges_of_verts_offsets[i + 1];
    unsigned nin = 0;
    for (unsigned j = f; j < e; ++j)
      if (crease_edges[edges_of_verts[j]])
        ++nin;
    corners[i] = (nin > 2);
  }
  return corners;
}

static unsigned* get_vert_class_dim(
    unsigned nverts,
    unsigned const* corner_verts,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_offsets,
    unsigned const* class_dim_of_edges)
{
  unsigned* class_dim = LOOP_MALLOC(unsigned, nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned f = edges_of_verts_offsets[i];
    unsigned e = edges_of_verts_offsets[i + 1];
    unsigned d = 3;
    for (unsigned j = f; j < e; ++j) {
      unsigned ed = class_dim_of_edges[edges_of_verts[j]];
      if (ed < d)
        d = ed;
    }
    class_dim[i] = (d - corner_verts[i]);
  }
  return class_dim;

}

void mesh_derive_class_dim(struct mesh* m, double crease_angle)
{
  unsigned dim = mesh_dim(m);
  unsigned* elem_class_dim = uints_filled(mesh_count(m, dim), dim);
  mesh_add_tag(m, dim, TAG_U32, "class_dim", 1, elem_class_dim);
  if (dim == 0)
    return;
  unsigned* boundary_sides = mesh_mark_part_boundary(m);
  unsigned nsides = mesh_count(m, dim - 1);
  unsigned* side_class_dim = LOOP_MALLOC(unsigned, nsides);
  for (unsigned i = 0; i < nsides; ++i)
    side_class_dim[i] = (dim - boundary_sides[i]);
  mesh_add_tag(m, dim - 1, TAG_U32, "class_dim", 1, side_class_dim);
  if (dim == 1) {
    loop_free(boundary_sides);
    return;
  }
  double* side_normals = get_boundary_side_normals(dim - 1, nsides,
      boundary_sides, mesh_ask_down(m, dim - 1, 0),
      mesh_find_tag(m, 0, "coordinates")->d.f64);
  unsigned* boundary_hinges = mesh_mark_down(m, dim - 1, dim - 2,
      boundary_sides);
  unsigned nhinges = mesh_count(m, dim - 2);
  double* hinge_angles = get_hinge_angles(nhinges,
      mesh_ask_up(m, dim - 2, dim - 1)->adj,
      mesh_ask_up(m, dim - 2, dim - 1)->offsets,
      boundary_sides, boundary_hinges, side_normals);
  loop_free(boundary_sides);
  loop_free(side_normals);
  unsigned* crease_hinges = mark_creases(nhinges, hinge_angles, crease_angle);
  loop_free(hinge_angles);
  unsigned* hinge_class_dim = LOOP_MALLOC(unsigned, nhinges);
  for (unsigned i = 0; i < nhinges; ++i)
    hinge_class_dim[i] = (dim - boundary_hinges[i] - crease_hinges[i]);
  loop_free(boundary_hinges);
  mesh_add_tag(m, dim - 2, TAG_U32, "class_dim", 1, hinge_class_dim);
  if (dim == 2) {
    loop_free(crease_hinges);
    return;
  }
  unsigned nverts = mesh_count(m, 0);
  unsigned* corner_verts = mark_corners(nverts,
      mesh_ask_up(m, 0, 1)->adj, mesh_ask_up(m, 0, 1)->offsets,
      crease_hinges);
  loop_free(crease_hinges);
  unsigned* vert_class_dim = get_vert_class_dim(nverts, corner_verts, 
      mesh_ask_up(m, 0, 1)->adj, mesh_ask_up(m, 0, 1)->offsets,
      hinge_class_dim);
  loop_free(corner_verts);
  mesh_add_tag(m, 0, TAG_U32, "class_dim", 1, vert_class_dim);
}

/* This algorithm amounts to finding connected components of a graph.
   It uses an iterative (non-recursive) depth-first search.
   It seems difficult to parallelize via "loop", and the use case
   is for serial pre-processing of under-specified meshes. */

static void set_equal_order_class_id(struct mesh* m, unsigned dim)
{
  assert(loop_size() == 1);
  unsigned n = mesh_count(m, dim);
  unsigned const* class_dim = mesh_find_tag(m, dim, "class_dim")->d.u32;
  unsigned const* down_class_dim = mesh_find_tag(m, dim - 1, "class_dim")->d.u32;
  unsigned const* down = mesh_ask_down(m, dim, dim - 1);
  unsigned degree = the_down_degrees[dim][dim - 1];
  unsigned const* up = mesh_ask_up(m, dim - 1, dim)->adj;
  unsigned const* up_offsets = mesh_ask_up(m, dim - 1, dim)->offsets;
  unsigned* class_id = uints_filled(n, INVALID);
  unsigned* stack = LOOP_MALLOC(unsigned, n);
  enum { WHITE, GRAY, BLACK };
  unsigned* state = uints_filled(n, WHITE);
  unsigned stack_n = 0;
  unsigned component = 0;
  for (unsigned i = 0; i < n; ++i)
    if (class_dim[i] != dim)
      state[i] = BLACK;
  for (unsigned i = 0; i < n; ++i) {
    if (state[i] != WHITE)
      continue;
    stack_n = 1;
    stack[0] = i;
    while (stack_n) {
      unsigned ent = stack[stack_n - 1];
      --stack_n;
      class_id[ent] = component;
      state[ent] = BLACK;
      for (unsigned j = 0; j < degree; ++j) {
        unsigned d = down[ent * degree + j];
        if (down_class_dim[d] != dim)
          continue;
        unsigned fu = up_offsets[d];
        unsigned eu = up_offsets[d + 1];
        for (unsigned k = fu; k < eu; ++k) {
          unsigned adj = up[k];
          if (state[adj] != WHITE)
            continue;
          state[adj] = GRAY;
          stack[stack_n] = adj;
          ++stack_n;
        }
      }
    }
    ++component;
  }
  loop_free(stack);
  loop_free(state);
  mesh_add_tag(m, dim, TAG_U32, "class_id", 1, class_id);
}

static void project_class_id_onto(struct mesh* m, unsigned low_dim)
{
  unsigned nlows = mesh_count(m, low_dim);
  unsigned const* high_class_dim = mesh_find_tag(m, low_dim + 1, "class_dim")->d.u32;
  unsigned const* high_class_id = mesh_find_tag(m, low_dim + 1, "class_id")->d.u32;
  unsigned const* low_class_dim = mesh_find_tag(m, low_dim, "class_dim")->d.u32;
  unsigned* low_class_id = mesh_find_tag(m, low_dim, "class_id")->d.u32;
  unsigned const* up = mesh_ask_up(m, low_dim, low_dim + 1)->adj;
  unsigned const* up_offsets = mesh_ask_up(m, low_dim, low_dim + 1)->offsets;
  for (unsigned i = 0; i < nlows; ++i) {
    unsigned f = up_offsets[i];
    unsigned e = up_offsets[i + 1];
    for (unsigned j = f; j < e; ++j) {
      unsigned high = up[j];
      if (high_class_dim[high] == low_class_dim[i])
        low_class_id[i] = high_class_id[high];
    }
  }
}

void mesh_derive_class_id(struct mesh* m)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned const* vert_class_dim = mesh_find_tag(m, 0, "class_dim")->d.u32;
  unsigned* vert_class_id = LOOP_MALLOC(unsigned, nverts);
  unsigned nmodel_verts = 0;
  for (unsigned i = 0; i < nverts; ++i)
    if (vert_class_dim[i] == 0)
      vert_class_id[i] = nmodel_verts++;
  mesh_add_tag(m, 0, TAG_U32, "class_id", 1, vert_class_id);
  unsigned dim = mesh_dim(m);
  for (unsigned d = 1; d <= dim; ++d)
    set_equal_order_class_id(m, d);
  for (unsigned dd = 1; dd <= dim; ++dd) {
    unsigned d = dim - dd;
    project_class_id_onto(m, d);
  }
}

void mesh_derive_model(struct mesh* m, double crease_angle)
{
  mesh_derive_class_dim(m, crease_angle);
  mesh_derive_class_id(m);
}
