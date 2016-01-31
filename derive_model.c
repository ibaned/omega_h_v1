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

LOOP_KERNEL(get_hinge_angle,
    unsigned const* sides_of_hinges,
    unsigned const* sides_of_hinges_offsets,
    unsigned const* boundary_sides,
    unsigned const* boundary_hinges,
    double const* side_normals,
    double* angles)

  if (!boundary_hinges[i]) {
    angles[i] = 0;
    return;
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

static double* get_hinge_angles(
    unsigned nhinges,
    unsigned const* sides_of_hinges,
    unsigned const* sides_of_hinges_offsets,
    unsigned const* boundary_sides,
    unsigned const* boundary_hinges,
    double const* side_normals)
{
  double* angles = LOOP_MALLOC(double, nhinges);
  LOOP_EXEC(get_hinge_angle, nhinges,
      sides_of_hinges, sides_of_hinges_offsets,
      boundary_sides, boundary_hinges, side_normals,
      angles);
  return angles;
}

LOOP_KERNEL(mark_crease,
    double const* hingle_angles,
    double crease_angle,
    unsigned* creases)
  creases[i] = hingle_angles[i] > crease_angle;
}

static unsigned* mark_creases(
    unsigned nhinges,
    double const* hingle_angles,
    double crease_angle)
{
  unsigned* creases = LOOP_MALLOC(unsigned, nhinges);
  LOOP_EXEC(mark_crease, nhinges, hingle_angles, crease_angle, creases);
  return creases;
}

LOOP_KERNEL(mark_corner,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_offsets,
    unsigned const* crease_edges,
    unsigned* corners)
  unsigned f = edges_of_verts_offsets[i];
  unsigned e = edges_of_verts_offsets[i + 1];
  unsigned nin = 0;
  for (unsigned j = f; j < e; ++j)
    if (crease_edges[edges_of_verts[j]])
      ++nin;
  corners[i] = (nin > 2);
}

static unsigned* mark_corners(
    unsigned nverts,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_offsets,
    unsigned const* crease_edges)
{
  unsigned* corners = LOOP_MALLOC(unsigned, nverts);
  LOOP_EXEC(mark_corner, nverts,
      edges_of_verts, edges_of_verts_offsets,
      crease_edges, corners);
  return corners;
}

LOOP_KERNEL(get_vert_class_dim_kern,
    unsigned const* corner_verts,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_offsets,
    unsigned const* class_dim_of_edges,
    unsigned* class_dim)
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

static unsigned* get_vert_class_dim(
    unsigned nverts,
    unsigned const* corner_verts,
    unsigned const* edges_of_verts,
    unsigned const* edges_of_verts_offsets,
    unsigned const* class_dim_of_edges)
{
  unsigned* class_dim = LOOP_MALLOC(unsigned, nverts);
  LOOP_EXEC(get_vert_class_dim_kern, nverts,
      corner_verts, edges_of_verts, edges_of_verts_offsets,
      class_dim_of_edges, class_dim);
  return class_dim;
}

LOOP_KERNEL(get_side_class_dim,
    unsigned dim,
    unsigned const* boundary_sides,
    unsigned* side_class_dim)
  side_class_dim[i] = (dim - boundary_sides[i]);
}

LOOP_KERNEL(get_hinge_class_dim,
    unsigned dim,
    unsigned const* boundary_hinges,
    unsigned const* crease_hinges,
    unsigned* hinge_class_dim)
  hinge_class_dim[i] = (dim - boundary_hinges[i] - crease_hinges[i]);
}

void mesh_derive_class_dim(struct mesh* m, double crease_angle)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  unsigned* elem_class_dim = uints_filled(nelems, dim);
  mesh_add_tag(m, dim, TAG_U32, "class_dim", 1, elem_class_dim);
  if (dim == 0)
    return;
  unsigned* boundary_sides = mesh_mark_part_boundary(m);
  unsigned nsides = mesh_count(m, dim - 1);
  unsigned* side_class_dim = LOOP_MALLOC(unsigned, nsides);
  LOOP_EXEC(get_side_class_dim, nsides, dim, boundary_sides, side_class_dim);
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
  LOOP_EXEC(get_hinge_class_dim, nhinges,
      dim, boundary_hinges, crease_hinges, hinge_class_dim);
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
   It seems difficult to parallelize via "loop", so we move data
   to the host and operate there. */

static void set_equal_order_class_id(struct mesh* m, unsigned dim)
{
  unsigned n = mesh_count(m, dim);
  unsigned ndown = mesh_count(m, dim - 1);
  unsigned* class_dim = LOOP_TO_HOST(unsigned,
      mesh_find_tag(m, dim, "class_dim")->d.u32, n);
  unsigned* down_class_dim = LOOP_TO_HOST(unsigned,
      mesh_find_tag(m, dim - 1, "class_dim")->d.u32, ndown);
  unsigned degree = the_down_degrees[dim][dim - 1];
  unsigned* down = LOOP_TO_HOST(unsigned,
      mesh_ask_down(m, dim, dim - 1), n * degree);
  unsigned* up = LOOP_TO_HOST(unsigned,
      mesh_ask_up(m, dim - 1, dim)->adj, n * degree);
  unsigned* up_offsets = LOOP_TO_HOST(unsigned,
      mesh_ask_up(m, dim - 1, dim)->offsets, ndown + 1);
  unsigned* class_id_dev = uints_filled(n, INVALID);
  unsigned* class_id = LOOP_TO_HOST(unsigned, class_id_dev, n);
  loop_free(class_id_dev);
  unsigned* stack = LOOP_HOST_MALLOC(unsigned, n);
  enum { WHITE, GRAY, BLACK };
  unsigned* state_dev = uints_filled(n, WHITE);
  unsigned* state = LOOP_TO_HOST(unsigned, state_dev, n);
  loop_free(state_dev);
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
  loop_host_free(class_dim);
  loop_host_free(down_class_dim);
  loop_host_free(down);
  loop_host_free(up);
  loop_host_free(up_offsets);
  loop_host_free(stack);
  loop_host_free(state);
  mesh_add_tag(m, dim, TAG_U32, "class_id", 1,
      LOOP_TO_DEVICE(unsigned, class_id, n));
  loop_host_free(class_id);
}

LOOP_KERNEL(project_class_kern,
    unsigned const* up_offsets,
    unsigned const* up,
    unsigned const* high_class_dim,
    unsigned const* high_class_id,
    unsigned const* low_class_dim,
    unsigned* low_class_id)
  unsigned f = up_offsets[i];
  unsigned e = up_offsets[i + 1];
  for (unsigned j = f; j < e; ++j) {
    unsigned high = up[j];
    if (high_class_dim[high] == low_class_dim[i])
      low_class_id[i] = high_class_id[high];
  }
}

static void project_class_id_onto(struct mesh* m, unsigned low_dim)
{
  unsigned nlows = mesh_count(m, low_dim);
  unsigned const* up_offsets = mesh_ask_up(m, low_dim, low_dim + 1)->offsets;
  unsigned const* up = mesh_ask_up(m, low_dim, low_dim + 1)->adj;
  unsigned const* high_class_dim =
    mesh_find_tag(m, low_dim + 1, "class_dim")->d.u32;
  unsigned const* high_class_id =
    mesh_find_tag(m, low_dim + 1, "class_id")->d.u32;
  unsigned const* low_class_dim =
    mesh_find_tag(m, low_dim, "class_dim")->d.u32;
  unsigned* low_class_id =
    mesh_find_tag(m, low_dim, "class_id")->d.u32;
  LOOP_EXEC(project_class_kern, nlows,
      up_offsets, up,
      high_class_dim, high_class_id,
      low_class_dim, low_class_id);
}

void mesh_derive_class_id(struct mesh* m)
{
  /* we'll do this part on the host as well, for laziness */
  {
    unsigned nverts = mesh_count(m, 0);
    unsigned* vert_class_dim = LOOP_TO_HOST(unsigned,
        mesh_find_tag(m, 0, "class_dim")->d.u32, nverts);
    unsigned* vert_class_id = LOOP_HOST_MALLOC(unsigned, nverts);
    unsigned nmodel_verts = 0;
    for (unsigned i = 0; i < nverts; ++i)
      if (vert_class_dim[i] == 0)
        vert_class_id[i] = nmodel_verts++;
    loop_host_free(vert_class_dim);
    mesh_add_tag(m, 0, TAG_U32, "class_id", 1,
        LOOP_TO_DEVICE(unsigned, vert_class_id, nverts));
    loop_host_free(vert_class_id);
  }
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
