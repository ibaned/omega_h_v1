#include "mesh.hpp"

#include <cassert>
#include <cstring>

#include "arrays.hpp"
#include "bridge_graph.hpp"
#include "derive_sides.hpp"
#include "dual.hpp"
#include "graph.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "parallel_mesh.hpp"
#include "reflect_down.hpp"
#include "star.hpp"
#include "tables.hpp"
#include "up_from_down.hpp"

struct up {
  unsigned* offsets;
  unsigned* adj;
  unsigned* directions;
};

struct mesh {
  unsigned elem_dim;
  enum mesh_rep rep;
  unsigned counts[4];
  unsigned* down[4][4];
  struct up* up[4][4];
  struct graph* star[4][4];
  unsigned* dual;
  struct tags tags[4];
  struct parallel_mesh* parallel;
};

static struct up* new_up(unsigned* offsets, unsigned* adj, unsigned* directions)
{
  struct up* u = LOOP_HOST_MALLOC(struct up, 1);
  u->offsets = offsets;
  u->adj = adj;
  u->directions = directions;
  return u;
}

static void free_up(struct up* u)
{
  if (!u)
    return;
  loop_free(u->offsets);
  loop_free(u->adj);
  loop_free(u->directions);
  loop_host_free(u);
}

struct mesh* new_mesh(unsigned elem_dim, enum mesh_rep rep, unsigned is_parallel)
{
  struct mesh* m = LOOP_HOST_MALLOC(struct mesh, 1);
  memset(m, 0, sizeof(*m));
  m->elem_dim = elem_dim;
  m->rep = rep;
  if (is_parallel)
    m->parallel = new_parallel_mesh();
  return m;
}

struct mesh* new_box_mesh(unsigned elem_dim)
{
  struct mesh* m = new_mesh(elem_dim, MESH_REDUCED, 0);
  unsigned nelems = the_box_nelems[elem_dim];
  unsigned nverts = the_box_nverts[elem_dim];
  mesh_set_ents(m, 0, nverts, 0);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* verts_of_elems = array_to_device(
      the_box_conns[elem_dim], verts_per_elem * nelems);
  mesh_set_ents(m, elem_dim, nelems, verts_of_elems);
  double* coords = array_to_device(
      the_box_coords[elem_dim], 3 * nverts);
  mesh_add_tag(m, 0, TAG_F64, "coordinates", 3, coords);
  return m;
}

unsigned mesh_dim(struct mesh* m)
{
  return m->elem_dim;
}

unsigned mesh_count(struct mesh* m, unsigned dim)
{
  if (dim > m->elem_dim)
    return 0;
  /* intermediate entities are counted when we derive their
   * vertices, so trigger this process if necessary
   */
  if (dim)
    mesh_ask_down(m, dim, 0);
  return m->counts[dim];
}

static void free_mesh_contents(struct mesh* m)
{
  for (unsigned high_dim = 0; high_dim <= m->elem_dim; ++high_dim)
    for (unsigned low_dim = 0; low_dim <= high_dim; ++low_dim) {
      loop_free(m->down[high_dim][low_dim]);
      free_up(m->up[low_dim][high_dim]);
      osh_free_graph(m->star[low_dim][high_dim]);
    }
  loop_free(m->dual);
  for (unsigned d = 0; d < 4; ++d)
    free_tags(&m->tags[d]);
  if (m->parallel)
    free_parallel_mesh(m->parallel);
}

void free_mesh(struct mesh* m)
{
  free_mesh_contents(m);
  loop_host_free(m);
}

struct const_tag* mesh_find_tag(struct mesh* m, unsigned dim,
    char const* name)
{
  return find_tag(&m->tags[dim], name);
}

static void set_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned* adj)
{
  assert( ! m->down[high_dim][low_dim]);
  m->down[high_dim][low_dim] = adj;
}

unsigned const* mesh_ask_down(struct mesh* m, unsigned high_dim, unsigned low_dim)
{
  assert(high_dim <= mesh_dim(m));
  assert(low_dim <= high_dim);
  if (m->down[high_dim][low_dim])
    return m->down[high_dim][low_dim];
  if (low_dim == high_dim) {
    /* waste memory to prevent algorithms from having to deal
       with equal-order cases separately */
    unsigned n = m->counts[high_dim];
    unsigned* lows_of_highs = uints_linear(n, 1);
    set_down(m, high_dim, low_dim, lows_of_highs);
  } else {
    if (low_dim > 0) {/* deriving intermediate downward adjacency */
      set_down(m, high_dim, low_dim,
          mesh_reflect_down(m, high_dim, low_dim));
    } else {/* deriving implicit entity to vertex connectivity */
      assert(mesh_get_rep(m) == MESH_REDUCED);
      if (high_dim == 1 && m->elem_dim == 3) { /* deriving edges in 3D */
        struct const_graph* verts_of_verts = mesh_ask_star(m, 0, m->elem_dim);
        unsigned nverts = m->counts[0];
        unsigned nedges;
        unsigned* verts_of_edges;
        bridge_graph(nverts, verts_of_verts->offsets, verts_of_verts->adj,
            &nedges, &verts_of_edges);
        mesh_set_ents(m, 1, nedges, verts_of_edges);
      } else { /* deriving sides (2D edges, 3D faces) */
        assert(high_dim == m->elem_dim - 1);
        unsigned const* elems_of_elems = mesh_ask_dual(m);
        unsigned nelems = m->counts[m->elem_dim];
        unsigned const* verts_of_elems = m->down[m->elem_dim][0];
        unsigned nsides;
        unsigned* elems_of_sides;
        unsigned* elem_side_of_sides;
        bridge_dual_graph(m->elem_dim, nelems, elems_of_elems,
            &nsides, &elems_of_sides, &elem_side_of_sides);
        unsigned* verts_of_sides = derive_sides(m->elem_dim, nsides,
            verts_of_elems, elems_of_sides, elem_side_of_sides);
        loop_free(elems_of_sides);
        loop_free(elem_side_of_sides);
        mesh_set_ents(m, high_dim, nsides, verts_of_sides);
      }
    }
  }
  return m->down[high_dim][low_dim];
}

static void set_up(struct mesh* m, unsigned low_dim, unsigned high_dim,
    struct up* adj)
{
  assert( ! m->up[low_dim][high_dim]);
  m->up[low_dim][high_dim] = adj;
}

struct const_up* mesh_ask_up(struct mesh* m, unsigned low_dim, unsigned high_dim)
{
  assert(low_dim <= high_dim);
  if (m->up[low_dim][high_dim])
    return reinterpret_cast<struct const_up*>(m->up[low_dim][high_dim]);
  if (low_dim == high_dim) {
    /* waste memory to prevent algorithms from having to deal
       with equal-order cases separately */
    unsigned n = m->counts[high_dim];
    unsigned* ones = uints_filled(n, 1);
    unsigned* offsets = uints_exscan(ones, n);
    loop_free(ones);
    unsigned* highs_of_lows = copy_array(offsets, n);
    unsigned* directions = uints_filled(n, 0);
    set_up(m, low_dim, high_dim, new_up(offsets, highs_of_lows, directions));
  } else {
    unsigned const* lows_of_highs = mesh_ask_down(m, high_dim, low_dim);
    unsigned nhighs = m->counts[high_dim];
    unsigned nlows = m->counts[low_dim];
    unsigned* offsets;
    unsigned* highs_of_lows;
    unsigned* directions;
    up_from_down(high_dim, low_dim, nhighs, nlows, lows_of_highs,
        &offsets, &highs_of_lows, &directions);
    set_up(m, low_dim, high_dim, new_up(offsets, highs_of_lows, directions));
  }
  return reinterpret_cast<struct const_up*>(m->up[low_dim][high_dim]);
}

static void set_star(struct mesh* m, unsigned low_dim, unsigned high_dim,
    struct graph* adj)
{
  assert( ! m->star[low_dim][high_dim]);
  m->star[low_dim][high_dim] = adj;
}

struct const_graph* mesh_ask_star(struct mesh* m, unsigned low_dim, unsigned high_dim)
{
  assert(low_dim <= high_dim);
  if (m->star[low_dim][high_dim])
    return reinterpret_cast<struct const_graph*>(m->star[low_dim][high_dim]);
  if (low_dim == high_dim) {
    /* waste memory to prevent algorithms from having to deal
       with equal-order cases separately */
    unsigned n = m->counts[low_dim];
    unsigned* offsets = uints_filled(n + 1, 0);
    set_star(m, low_dim, high_dim, osh_new_graph(offsets, 0));
  } else {
    unsigned* offsets;
    unsigned* adj;
    mesh_get_star(m, low_dim, high_dim, &offsets, &adj);
    set_star(m, low_dim, high_dim, osh_new_graph(offsets, adj));
  }
  return reinterpret_cast<struct const_graph*>(m->star[low_dim][high_dim]);
}

static void set_dual(struct mesh* m, unsigned* adj)
{
  assert( ! m->dual);
  m->dual = adj;
}

unsigned const* mesh_ask_dual(struct mesh* m)
{
  if (!m->dual) {
    if (mesh_has_dim(m, mesh_dim(m) - 1))
      set_dual(m, mesh_get_dual_from_sides(m));
    else
      set_dual(m, mesh_get_dual_from_verts(m));
  }
  return m->dual;
}

void mesh_set_ents(struct mesh* m, unsigned dim, unsigned n, unsigned* verts)
{
  m->counts[dim] = n;
  set_down(m, dim, 0, verts);
}

struct const_tag* mesh_add_tag(struct mesh* m, unsigned dim, enum tag_type type,
    char const* name, unsigned ncomps, void* data)
{
  return add_tag(&m->tags[dim], type, name, ncomps, data);
}

void mesh_free_tag(struct mesh* m, unsigned dim, char const* name)
{
  remove_tag(&m->tags[dim], name);
}

unsigned mesh_count_tags(struct mesh* m, unsigned dim)
{
  return count_tags(&m->tags[dim]);
}

struct const_tag* mesh_get_tag(struct mesh* m, unsigned dim, unsigned i)
{
  return get_tag(&m->tags[dim], i);
}

struct tags* mesh_tags(struct mesh* m, unsigned dim)
{
  return &m->tags[dim];
}

unsigned mesh_has_dim(struct mesh* m, unsigned dim)
{
  return dim == 0 || m->down[dim][0] != 0;
}

enum mesh_rep mesh_get_rep(struct mesh* m)
{
  return m->rep;
}

void mesh_set_rep(struct mesh* m, enum mesh_rep rep)
{
  m->rep = rep;
}

unsigned mesh_is_parallel(struct mesh* m)
{
  return m->parallel != 0;
}

struct parallel_mesh* mesh_parallel(struct mesh* m)
{
  return m->parallel;
}

void mesh_make_parallel(struct mesh* m)
{
  assert(!mesh_is_parallel(m));
  m->parallel = new_parallel_mesh();
  for (unsigned d = 0; d <= mesh_dim(m); ++d)
    if (mesh_has_dim(m, d))
      mesh_set_globals(m, d, ulongs_linear(mesh_count(m, d), 1));
}

void overwrite_mesh(struct mesh* old, struct mesh* with)
{
  free_mesh_contents(old);
  *old = *with;
  loop_host_free(with);
}
