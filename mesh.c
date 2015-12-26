#include "mesh.h"

#include <assert.h>
#include <string.h>

#include "arrays.h"
#include "bridge_graph.h"
#include "derive_sides.h"
#include "graph.h"
#include "ints.h"
#include "loop.h"
#include "parallel_mesh.h"
#include "reflect_down.h"
#include "star.h"
#include "tables.h"
#include "up_from_down.h"

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

struct mesh* new_mesh(unsigned elem_dim)
{
  struct mesh* m = LOOP_HOST_MALLOC(struct mesh, 1);
  memset(m, 0, sizeof(*m));
  m->elem_dim = elem_dim;
  m->parallel = new_parallel_mesh(m);
  return m;
}

struct mesh* new_box_mesh(unsigned elem_dim)
{
  struct mesh* m = new_mesh(elem_dim);
  unsigned nelems = the_box_nelems[elem_dim];
  unsigned nverts = the_box_nverts[elem_dim];
  mesh_set_ents(m, 0, nverts, 0);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned nbytes = sizeof(unsigned) * verts_per_elem * nelems;
  unsigned* verts_of_elems = LOOP_MALLOC(unsigned, verts_per_elem * nelems);
  memcpy(verts_of_elems, the_box_conns[elem_dim], nbytes);
  double* coords = LOOP_MALLOC(double, 3 * nverts);
  nbytes = sizeof(double) * 3 * nverts;
  memcpy(coords, the_box_coords[elem_dim], nbytes);
  mesh_set_ents(m, elem_dim, nelems, verts_of_elems);
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

void free_mesh(struct mesh* m)
{
  for (unsigned high_dim = 0; high_dim <= m->elem_dim; ++high_dim)
    for (unsigned low_dim = 0; low_dim <= high_dim; ++low_dim) {
      loop_free(m->down[high_dim][low_dim]);
      free_up(m->up[low_dim][high_dim]);
      free_graph(m->star[low_dim][high_dim]);
    }
  loop_free(m->dual);
  for (unsigned d = 0; d < 4; ++d)
    free_tags(&m->tags[d]);
  free_parallel_mesh(m->parallel);
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
    if (low_dim) {/* deriving intermediate downward adjacency */
      unsigned nhighs = m->counts[high_dim];
      unsigned const* verts_of_highs = mesh_ask_down(m, high_dim, 0);
      struct const_up* lows_of_verts = mesh_ask_up(m, 0, low_dim);
      unsigned* lows_of_highs = reflect_down(high_dim, low_dim, nhighs,
          verts_of_highs, lows_of_verts->offsets, lows_of_verts->adj);
      set_down(m, high_dim, low_dim, lows_of_highs);
    } else {/* deriving intermediate entities (entity to vertex connectivity) */
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
    return (struct const_up*) m->up[low_dim][high_dim];
  if (low_dim == high_dim) {
    /* waste memory to prevent algorithms from having to deal
       with equal-order cases separately */
    unsigned n = m->counts[high_dim];
    unsigned* ones = uints_filled(n, 1);
    unsigned* offsets = uints_exscan(ones, n);
    loop_free(ones);
    unsigned* highs_of_lows = uints_copy(offsets, n);
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
  return (struct const_up*) m->up[low_dim][high_dim];
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
    return (struct const_graph*) m->star[low_dim][high_dim];
  if (low_dim == high_dim) {
    /* waste memory to prevent algorithms from having to deal
       with equal-order cases separately */
    unsigned n = m->counts[low_dim];
    unsigned* offsets = uints_filled(n + 1, 0);
    set_star(m, low_dim, high_dim, new_graph(offsets, 0));
  } else {
    unsigned const* lows_of_highs = mesh_ask_down(m, high_dim, low_dim);
    struct const_up* highs_of_lows = mesh_ask_up(m, low_dim, high_dim);
    unsigned nlows = m->counts[low_dim];
    unsigned* offsets;
    unsigned* adj;
    get_star(low_dim, high_dim, nlows, highs_of_lows->offsets, highs_of_lows->adj,
        lows_of_highs, &offsets, &adj);
    set_star(m, low_dim, high_dim, new_graph(offsets, adj));
  }
  return (struct const_graph*) m->star[low_dim][high_dim];
}

static void set_dual(struct mesh* m, unsigned* adj)
{
  assert( ! m->dual);
  m->dual = adj;
}

unsigned const* mesh_ask_dual(struct mesh* m)
{
  if (m->dual)
    return m->dual;
  unsigned nelems = m->counts[m->elem_dim];
  unsigned const* verts_of_elems = m->down[m->elem_dim][0];
  struct const_up* elems_of_verts = mesh_ask_up(m, 0, m->elem_dim);
  unsigned* elems_of_elems = get_dual(m->elem_dim, nelems, verts_of_elems,
      elems_of_verts->offsets, elems_of_verts->adj);
  set_dual(m, elems_of_elems);
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

struct parallel_mesh* mesh_parallel(struct mesh* m)
{
  return m->parallel;
}

enum mesh_rep mesh_get_rep(struct mesh* m)
{
  return m->rep;
}

void mesh_set_rep(struct mesh* m, enum mesh_rep r)
{
  m->rep = r;
}

unsigned mesh_is_parallel(struct mesh* m)
{
  return m->parallel != 0;
}

void mesh_set_parallel(struct mesh* m, unsigned yn)
{
  if (yn && !m->parallel)
    m->parallel = new_parallel_mesh(m);
  else if (!yn && !m->parallel) {
    free_parallel_mesh(m->parallel);
    m->parallel = 0;
  }
}
