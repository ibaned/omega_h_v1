#include "ghost_mesh.h"

#include <assert.h>
#include <string.h>

#include "arrays.h"
#include "close_partition.h"
#include "exchanger.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "migrate_mesh.h"
#include "parallel_mesh.h"
#include "subset.h"

/* this function is analogous to
   get_vert_use_owners_of_elems
   in close_partition.c.
   each owned vertex in the old mesh is
   given a set of adjacent elements.
   that set is the correct global set,
   accounting for all partitions.
   each element is, as usual, identified
   by its owner.
   the main difference in this code is
   that adjacent elements have to be
   collected from all the partitions. */

static void get_elem_use_owners_of_verts(
    struct mesh* m,
    unsigned** p_use_own_ranks,
    unsigned** p_use_own_ids,
    unsigned** p_uses_of_verts_offsets)
{
  unsigned dim = mesh_dim(m);
  unsigned const* elems_of_verts_offsets =
    mesh_ask_up(m, 0, dim)->offsets;
  unsigned const* elems_of_verts =
    mesh_ask_up(m, 0, dim)->adj;
  unsigned const* elem_own_ranks = mesh_ask_own_ranks(m, dim);
  unsigned const* elem_own_ids = mesh_ask_own_ids(m, dim);
  unsigned nverts = mesh_count(m, 0);
  unsigned nuses_in = elems_of_verts_offsets[nverts];
  unsigned* use_ranks_in = LOOP_MALLOC(unsigned, nuses_in);
  unsigned* use_ids_in = LOOP_MALLOC(unsigned, nuses_in);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned f = elems_of_verts_offsets[i];
    unsigned e = elems_of_verts_offsets[i + 1];
    for (unsigned j = f; j < e; ++j) {
      unsigned elem = elems_of_verts[j];
      use_ranks_in[j] = elem_own_ranks[elem];
      use_ids_in[j] = elem_own_ids[elem];
    }
  }
  unsigned const* vert_own_ranks = mesh_ask_own_ranks(m, 0);
  unsigned const* vert_own_ids = mesh_ask_own_ids(m, 0);
  unsigned* dest_ranks = uints_expand(nverts, 1,
      vert_own_ranks, elems_of_verts_offsets);
  unsigned* dest_ids = uints_expand(nverts, 1,
      vert_own_ids, elems_of_verts_offsets);
  struct exchanger* ex = new_exchanger(nuses_in, dest_ranks);
  loop_free(dest_ranks);
  set_exchanger_dests(ex, nverts, dest_ids);
  loop_free(dest_ids);
  *p_use_own_ranks = exchange_uints(ex, 1, use_ranks_in,
      EX_FOR, EX_ITEM);
  loop_free(use_ranks_in);
  *p_use_own_ids = exchange_uints(ex, 1, use_ids_in,
      EX_FOR, EX_ITEM);
  loop_free(use_ids_in);
  *p_uses_of_verts_offsets = uints_copy(
      ex->items_of_roots_offsets[EX_REV], nverts + 1);
  free_exchanger(ex);
}

/* the set of entities on this partition in the
   new mesh, identified by their owners in the old mesh */

struct resident {
  unsigned n;
  unsigned padding__;
  unsigned* ranks;
  unsigned* ids;
};

/* the set of entities B adjacent to a given set of entities A,
   identified by the owners of B in the old mesh and
   offsets from the A group to the B group */

struct uses {
  unsigned* offsets;
  unsigned* ranks;
  unsigned* ids;
};

enum ghost_type { VERT, ELEM };

struct ghost_state {
  unsigned nown[2];
  /* the entities currently resident on this part */
  struct resident resident[2];
  /* the "opposite" entities adjacent to the resident entities */
  struct uses res_uses[2];
  /* the opposite entities adjacent to the old owners
     (used by push_use_owners to set up res_uses) */
  struct uses own_uses[2];
};

static unsigned ghost_dim(struct mesh* m, enum ghost_type t)
{
  switch (t) {
    case VERT: return 0;
    case ELEM: return mesh_dim(m);
    default: return 42;
  }
}

static enum ghost_type ghost_opp(enum ghost_type t)
{
  switch (t) {
    case VERT: return ELEM;
    case ELEM: return VERT;
    default: return VERT;
  }
}

static void init_ghost_resident(struct ghost_state* s,
    struct mesh* m, enum ghost_type t)
{
  unsigned d = ghost_dim(m, t);
  unsigned n = mesh_count(m, d);
  s->resident[t].n = n;
  s->resident[t].ranks = uints_copy(mesh_ask_own_ranks(m, d), n);
  s->resident[t].ids = uints_copy(mesh_ask_own_ids(m, d), n);
}

static void init_ghost_uses(struct ghost_state* s,
    struct mesh* m, enum ghost_type t)
{
  switch (t) {
    case VERT: get_elem_use_owners_of_verts(m,
          &s->own_uses[t].ranks,
          &s->own_uses[t].ids,
          &s->own_uses[t].offsets);
      break;
    case ELEM: get_down_use_owners(m,
          mesh_dim(m), 0,
          &s->own_uses[t].ranks,
          &s->own_uses[t].ids,
          &s->own_uses[t].offsets);
      break;
  }
}

static void free_resident(struct resident* r)
{
  loop_free(r->ranks);
  loop_free(r->ids);
}

static void free_uses(struct uses* u)
{
  loop_free(u->ranks);
  loop_free(u->ids);
  loop_free(u->offsets);
}

static void init_ghosts(struct ghost_state* s, struct mesh* m)
{
  s->nown[VERT] = mesh_count(m, 0);
  s->nown[ELEM] = mesh_count(m, mesh_dim(m));
  /* the old mesh already tells us which vertices are resident
     in a non-ghosted mesh, we'll start from there */
  init_ghost_resident(s, m, VERT);
  init_ghost_uses(s, m, VERT);
  init_ghost_uses(s, m, ELEM);
}

/* figure out the adjacencies of resident entities based
   on the adjacencies of their owners */

static void push_ghosts(struct ghost_state* s, enum ghost_type t)
{
  struct exchanger* push = make_reverse_exchanger(s->nown[t],
      s->resident[t].n, s->resident[t].ranks, s->resident[t].ids);
  free_uses(&s->res_uses[t]);
  push_use_owners(push,
      s->own_uses[t].ranks, s->own_uses[t].ids, s->own_uses[t].offsets,
      &s->res_uses[t].ranks, &s->res_uses[t].ids, &s->res_uses[t].offsets);
  free_exchanger(push);
}

/* figure out the resident entities of type B
   based on the resident entities of type A
   and their adjacent type B entities */

static void close_ghosts(struct ghost_state* s, enum ghost_type t)
{
  enum ghost_type ot = ghost_opp(t);
  free_resident(&s->resident[ot]);
  close_partition(s->resident[t].n, s->nown[ot],
      s->res_uses[t].offsets, s->res_uses[t].ranks, s->res_uses[t].ids,
      &s->resident[ot].n, &s->resident[ot].ranks, &s->resident[ot].ids);
}

static void free_ghosts(struct ghost_state* s)
{
  for (unsigned i = 0; i < 2; ++i) {
    free_uses(&s->own_uses[i]);
    free_uses(&s->res_uses[i]);
  }
  free_resident(&s->resident[VERT]);
  /* we don't free the resident elements,
     that gets input to mesh_migrate */
}

void ghost_mesh(struct mesh** p_m, unsigned nlayers)
{
  assert(mesh_ghost_layers(*p_m) == 0);
  if (nlayers == 0)
    return;
  struct ghost_state s;
  memset(&s, 0, sizeof(s));
  init_ghosts(&s, *p_m);
  push_ghosts(&s, VERT);
  close_ghosts(&s, VERT);
  /* now we have the resident elements for 1 layer of ghosting */
  for (unsigned i = 1; i < nlayers; ++i) {
    push_ghosts(&s, ELEM);
    close_ghosts(&s, ELEM);
    /* now we have the resident vertices for (i) layers of ghosting */
    push_ghosts(&s, VERT);
    close_ghosts(&s, VERT);
    /* now we have the resident vertices for (i + 1) layers of ghosting */
  }
  free_ghosts(&s); /* deletes all but resident elements */
/* if we let the owners be generated based on global numbers,
   then the ghost layers would not be guaranteed to be un-owned.
   since we are using element ownership to identify ghost layers,
   we need to strictly preserve the same owner rank for entities
   as we increase their copies through ghosting */
  for (unsigned d = 0; d <= mesh_dim(*p_m); ++d)
    if (mesh_has_dim(*p_m, d))
      mesh_tag_own_rank(*p_m, d);
  migrate_mesh(p_m, s.resident[ELEM].n,
      s.resident[ELEM].ranks, s.resident[ELEM].ids);
  free_resident(&s.resident[ELEM]);
  mesh_set_ghost_layers(*p_m, nlayers);
}

void unghost_mesh(struct mesh** p_m)
{
  struct mesh* m = *p_m;
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  unsigned* owned_elems = mesh_get_owned(m, dim);
  unsigned* offsets = uints_exscan(owned_elems, nelems);
  loop_free(owned_elems);
  struct mesh* ugm = subset_mesh(m, dim, offsets);
  loop_free(offsets);
  free_mesh(m);
  *p_m = ugm;
}

void mesh_ensure_ghosting(struct mesh** p_m, unsigned nlayers)
{
  if (nlayers == mesh_ghost_layers(*p_m))
    return;
  if (mesh_ghost_layers(*p_m))
    unghost_mesh(p_m);
  if (nlayers)
    ghost_mesh(p_m, nlayers);
}
