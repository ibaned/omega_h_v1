#include "derive_model.hpp"

#include <cassert>
#include <cmath>

#include "algebra.hpp"
#include "arrays.hpp"
#include "bfs.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mark.hpp"
#include "mesh.hpp"
#include "tables.hpp"

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

LOOP_KERNEL(count_boundary_graph,
    unsigned const* eq_offsets,
    unsigned const* bridges_of_ents,
    unsigned bridges_per_ent,
    unsigned const* bridges,
    unsigned* degrees)
  if (eq_offsets[i] == eq_offsets[i + 1])
    return;
  unsigned d = 0;
  unsigned const* bridges_of_ent = bridges_of_ents + i * bridges_per_ent;
  for (unsigned j = 0; j < bridges_per_ent; ++j)
    if (bridges[bridges_of_ent[j]])
      ++d;
  degrees[eq_offsets[i]] = d;
}

LOOP_KERNEL(fill_boundary_graph,
    unsigned const* eq_offsets,
    unsigned const* bridges_of_ents,
    unsigned bridges_per_ent,
    unsigned const* bridges,
    unsigned const* ents_of_bridges_offsets,
    unsigned const* ents_of_bridges,
    unsigned const* offsets,
    unsigned* adj)
  if (eq_offsets[i] == eq_offsets[i + 1])
    return;
  unsigned eq = eq_offsets[i];
  unsigned k = offsets[eq];
  unsigned const* bridges_of_ent = bridges_of_ents + i * bridges_per_ent;
  for (unsigned j = 0; j < bridges_per_ent; ++j) {
    unsigned bridge = bridges_of_ent[j];
    if (bridges[bridges_of_ent[j]]) {
      unsigned a = ents_of_bridges_offsets[bridge];
      unsigned b = ents_of_bridges_offsets[bridge + 1];
      unsigned l;
      for (l = a; l < b; ++l) {
        unsigned other = ents_of_bridges[l];
        if (other == i)
          continue;
        if (eq_offsets[other] != eq_offsets[other + 1]) {
          adj[k++] = eq_offsets[other];
          break;
        }
      }
      assert(l < b);
    }
  }
}

static void form_boundary_graph(struct mesh* m, unsigned dim,
    unsigned** p_eq_offsets, unsigned** p_offsets, unsigned** p_adj)
{
  assert(dim > 0);
  unsigned nents = mesh_count(m, dim);
  unsigned* eq_ents = mesh_mark_class(m, dim, dim, INVALID);
  unsigned* eq_offsets = uints_exscan(eq_ents, nents);
  unsigned neqs = uints_at(eq_offsets, nents);
  loop_free(eq_ents);
  unsigned* bridges = mesh_mark_class(m, dim - 1, dim, INVALID);
  unsigned bridges_per_ent = the_down_degrees[dim][dim - 1];
  unsigned const* bridges_of_ents = mesh_ask_down(m, dim, dim - 1);
  unsigned* degrees = LOOP_MALLOC(unsigned, neqs);
  LOOP_EXEC(count_boundary_graph, nents, eq_offsets, bridges_of_ents,
      bridges_per_ent, bridges, degrees);
  unsigned* offsets = uints_exscan(degrees, neqs);
  loop_free(degrees);
  unsigned nadj = uints_at(offsets, neqs);
  unsigned* adj = LOOP_MALLOC(unsigned, nadj);
  unsigned const* ents_of_bridges_offsets =
    mesh_ask_up(m, dim - 1, dim)->offsets;
  unsigned const* ents_of_bridges =
    mesh_ask_up(m, dim - 1, dim)->adj;
  LOOP_EXEC(fill_boundary_graph, nents, eq_offsets, bridges_of_ents,
      bridges_per_ent, bridges, ents_of_bridges_offsets, ents_of_bridges,
      offsets, adj);
  loop_free(bridges);
  *p_eq_offsets = eq_offsets;
  *p_offsets = array_to_host(offsets, neqs + 1);
  loop_free(offsets);
  *p_adj = array_to_host(adj, nadj);
  loop_free(adj);
}

LOOP_KERNEL(extract_eq_class_id,
    unsigned const* eq_offsets,
    unsigned const* comp,
    unsigned* class_id)
  if (eq_offsets[i] != eq_offsets[i + 1])
    class_id[i] = comp[eq_offsets[i]];
}

static void set_equal_order_class_id(struct mesh* m, unsigned dim)
{
  unsigned* eq_offsets;
  unsigned* offsets;
  unsigned* adj;
  form_boundary_graph(m, dim, &eq_offsets, &offsets, &adj);
  unsigned nents = mesh_count(m, dim);
  unsigned neqs = uints_at(eq_offsets, nents);
  unsigned* host_comp = LOOP_HOST_MALLOC(unsigned, neqs);
  connected_components(neqs, offsets, adj, host_comp);
  loop_host_free(offsets);
  loop_host_free(adj);
  unsigned* comp = array_to_device(host_comp, neqs);
  loop_host_free(host_comp);
  unsigned* class_id = LOOP_MALLOC(unsigned, nents);
  LOOP_EXEC(extract_eq_class_id, nents, eq_offsets, comp, class_id);
  loop_free(comp);
  loop_free(eq_offsets);
  mesh_add_tag(m, dim, TAG_U32, "class_id", 1, class_id);
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
    unsigned* vert_class_dim = array_to_host(
        mesh_find_tag(m, 0, "class_dim")->d.u32, nverts);
    unsigned* vert_class_id = LOOP_HOST_MALLOC(unsigned, nverts);
    unsigned nmodel_verts = 0;
    for (unsigned i = 0; i < nverts; ++i)
      if (vert_class_dim[i] == 0)
        vert_class_id[i] = nmodel_verts++;
    loop_host_free(vert_class_dim);
    mesh_add_tag(m, 0, TAG_U32, "class_id", 1,
        array_to_device(vert_class_id, nverts));
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
