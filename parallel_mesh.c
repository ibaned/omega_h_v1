#include "parallel_mesh.h"

#include <assert.h>
#include <string.h>

#include "algebra.h"
#include "comm.h"
#include "doubles.h"
#include "exchanger.h"
#include "global.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "owners_from_global.h"

struct parallel_mesh {
  struct mesh* m;
  unsigned* own_ranks[4];
  unsigned* own_ids[4];
  struct exchanger* ex[4];
};

struct parallel_mesh* new_parallel_mesh(struct mesh* m)
{
  struct parallel_mesh* pm = LOOP_HOST_MALLOC(struct parallel_mesh, 1);
  memset(pm, 0, sizeof(*pm));
  pm->m = m;
  return pm;
}

void free_parallel_mesh(struct parallel_mesh* pm)
{
  for (unsigned i = 0; i < 4; ++i) {
    loop_free(pm->own_ranks[i]);
    loop_free(pm->own_ids[i]);
    if (pm->ex[i])
      free_exchanger(pm->ex[i]);
  }
  loop_host_free(pm);
}

unsigned long const* mesh_ask_global(struct mesh* m, unsigned dim)
{
  return mesh_find_tag(m, dim, "global_number")->d.u64;
}

static void ask_owners(struct mesh* m, unsigned dim)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  if (pm->own_ranks[dim])
    return;
  assert(mesh_find_tag(m, dim, "global_number"));
  owners_from_global(mesh_count(m, dim),
      mesh_find_tag(m, dim, "global_number")->d.u64,
      &pm->own_ranks[dim], &pm->own_ids[dim]);
  return;
}

unsigned const* mesh_ask_own_ranks(struct mesh* m, unsigned dim)
{
  ask_owners(m, dim);
  return mesh_parallel(m)->own_ranks[dim];
}

unsigned const* mesh_ask_own_ids(struct mesh* m, unsigned dim)
{
  ask_owners(m, dim);
  return mesh_parallel(m)->own_ids[dim];
}

struct exchanger* mesh_ask_exchanger(struct mesh* m, unsigned dim)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  if (!pm->ex[dim]) {
    unsigned const* own_ranks = mesh_ask_own_ranks(m, dim);
    unsigned const* own_ids = mesh_ask_own_ids(m, dim);
    unsigned n = mesh_count(m, dim);
    pm->ex[dim] = new_exchanger(n, own_ranks);
    set_exchanger_dests(pm->ex[dim], n, own_ids);
  }
  return pm->ex[dim];
}

static unsigned long* global_from_owners(
    struct exchanger* ex,
    unsigned const* own_ranks)
{
  unsigned n = ex->nitems[EX_REV];
  unsigned* owned = LOOP_MALLOC(unsigned, n);
  unsigned rank = comm_rank();
  for (unsigned i = 0; i < n; ++i)
    owned[i] = (own_ranks[i] == rank);
  unsigned* offsets = uints_exscan(owned, n);
  unsigned long* local_globals = globalize_offsets(offsets, n);
  unsigned long* globals = exchange_ulongs(ex, 1, local_globals,
      EX_REV, EX_ROOT);
  loop_free(local_globals);
  return globals;
}

void mesh_global_renumber(struct mesh* m, unsigned dim)
{
  unsigned long* new_global = global_from_owners(
      mesh_ask_exchanger(m, dim),
      mesh_ask_own_ranks(m, dim));
  if (mesh_find_tag(m, dim, "global_number"))
    mesh_free_tag(m, dim, "global_number");
  mesh_add_tag(m, dim, TAG_U64, "global_number", 1, new_global);
}

void mesh_conform_tag(struct mesh* m, unsigned dim, const char* name)
{
  struct const_tag* t = mesh_find_tag(m, dim, name);
  struct exchanger* ex = mesh_ask_exchanger(m, dim);
  exchange_tag(ex, t, mesh_tags(m, dim), EX_REV, EX_ROOT);
}

void mesh_accumulate_tag(struct mesh* m, unsigned dim, const char* name)
{
  struct const_tag* t = mesh_find_tag(m, dim, name);
  assert(t->type == TAG_F64);
  struct exchanger* ex = mesh_ask_exchanger(m, dim);
  double* in = exchange_doubles(ex, t->ncomps, t->d.f64, EX_FOR, EX_ITEM);
  unsigned nowners = mesh_count(m, dim);
  unsigned const* copies_of_owners_offsets =
    ex->items_of_roots_offsets[EX_REV];
  double* out = doubles_filled(nowners * t->ncomps, 0.0);
  for (unsigned i = 0; i < nowners; ++i) {
    unsigned f = copies_of_owners_offsets[i];
    unsigned e = copies_of_owners_offsets[i + 1];
    for (unsigned j = f; j < e; ++j)
      add_vectors(
          in + j * t->ncomps,
          out + i * t->ncomps,
          out + i * t->ncomps,
          t->ncomps);
  }
  loop_free(in);
  modify_tag(mesh_tags(m, dim), t->name, out);
}
