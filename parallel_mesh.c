#include "parallel_mesh.h"

#include <assert.h>
#include <string.h>

#include "exchanger.h"
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
  return mesh_parallel(m)->own_ranks[dim];
}

struct exchanger* mesh_ask_exchanger(struct mesh* m, unsigned dim)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  if (!pm->ex[dim]) {
    unsigned const* own_ranks = mesh_ask_own_ranks(m, dim);
    unsigned const* own_ids = mesh_ask_own_ids(m, dim);
    unsigned n = mesh_count(m, dim);
    pm->ex[dim] = new_exchanger(n, n, own_ranks, own_ids);
  }
  return pm->ex[dim];
}
