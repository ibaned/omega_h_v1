#include "parallel_mesh.hpp"

#include <cassert>
#include <cstring>

#include "arrays.hpp"
#include "algebra.hpp"
#include "bcast.hpp"
#include "comm.hpp"
#include "doubles.hpp"
#include "exchanger.hpp"
#include "global.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "owners_from_global.hpp"
#include "parallel_inertial_bisect.hpp"
#include "vtk_io.hpp"

namespace omega_h {

struct parallel_mesh {
  unsigned long* globals[4];
  unsigned* own_ranks[4];
  unsigned* own_ids[4];
  struct exchanger* ex[4];
  unsigned nghost_layers;
};

struct parallel_mesh* new_parallel_mesh(void)
{
  struct parallel_mesh* pm = LOOP_HOST_MALLOC(struct parallel_mesh, 1);
  memset(pm, 0, sizeof(*pm));
  return pm;
}

static void invalidate_exchanger(struct parallel_mesh* pm, unsigned dim)
{
  if (pm->ex[dim])
    free_exchanger(pm->ex[dim]);
  pm->ex[dim] = 0;
}

static void invalidate_ids(struct parallel_mesh* pm, unsigned dim)
{
  loop_free(pm->own_ids[dim]);
  pm->own_ids[dim] = 0;
  invalidate_exchanger(pm, dim);
}

static void invalidate_ranks(struct parallel_mesh* pm, unsigned dim)
{
  loop_free(pm->own_ranks[dim]);
  pm->own_ranks[dim] = 0;
  invalidate_ids(pm, dim);
}

static void invalidate_globals(struct parallel_mesh* pm, unsigned dim)
{
  loop_free(pm->globals[dim]);
  pm->globals[dim] = 0;
  invalidate_ranks(pm, dim);
}

void free_parallel_mesh(struct parallel_mesh* pm)
{
  for (unsigned i = 0; i < 4; ++i)
    invalidate_globals(pm, i);
  loop_host_free(pm);
}

unsigned long const* mesh_ask_globals(struct mesh* m, unsigned dim)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  return pm->globals[dim];
}

unsigned const* mesh_ask_own_ranks(struct mesh* m, unsigned dim)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  if (!pm->own_ranks[dim])
    owners_from_global(mesh_count(m, dim), mesh_ask_globals(m, dim),
        &pm->own_ranks[dim], &pm->own_ids[dim]);
  return pm->own_ranks[dim];
}

unsigned const* mesh_ask_own_ids(struct mesh* m, unsigned dim)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  if (!pm->own_ids[dim])
    own_idxs_from_global(mesh_count(m, dim), mesh_ask_globals(m, dim),
        pm->own_ranks[dim], &pm->own_ids[dim]);
  return pm->own_ids[dim];
}

struct exchanger* mesh_ask_exchanger(struct mesh* m, unsigned dim)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  if (!pm->ex[dim]) {
    unsigned const* own_ranks = mesh_ask_own_ranks(m, dim);
    unsigned const* own_ids = mesh_ask_own_ids(m, dim);
    unsigned n = mesh_count(m, dim);
    pm->ex[dim] = make_reverse_exchanger(n, n, own_ranks, own_ids);
  }
  return pm->ex[dim];
}

static unsigned long* global_from_owners(
    struct exchanger* ex,
    unsigned const* owned)
{
  unsigned nowners = ex->nroots[EX_FOR];
  unsigned* offsets = uints_exscan(owned, nowners);
  unsigned long* local_globals = globalize_offsets(offsets, nowners);
  loop_free(offsets);
  unsigned long* globals = exchange(ex, 1, local_globals,
      EX_FOR, EX_ROOT);
  loop_free(local_globals);
  return globals;
}

void mesh_global_renumber(struct mesh* m, unsigned dim)
{
  unsigned* owned = mesh_get_owned(m, dim);
  unsigned long* new_globals = global_from_owners(
      mesh_ask_exchanger(m, dim), owned);
  loop_free(owned);
  struct parallel_mesh* pm = mesh_parallel(m);
  if (pm->globals[dim])
    loop_free(pm->globals[dim]);
  pm->globals[dim] = new_globals;
}

template <typename T>
void mesh_conform_array(struct mesh* m, unsigned dim, unsigned width,
    T** a)
{
  if (!mesh_is_parallel(m))
    return;
  T* in = *a;
  T* out = exchange<T>(mesh_ask_exchanger(m, dim), width, in,
      EX_FOR, EX_ROOT);
  loop_free(in);
  *a = out;
}

template void mesh_conform_array(struct mesh* m, unsigned dim, unsigned width,
    double** a);
template void mesh_conform_array(struct mesh* m, unsigned dim, unsigned width,
    unsigned** a);
template void mesh_conform_array(struct mesh* m, unsigned dim, unsigned width,
    unsigned long** a);

void mesh_conform_tag(struct mesh* m, unsigned dim, const char* name)
{
  if (!mesh_is_parallel(m))
    return;
  struct const_tag* t = mesh_find_tag(m, dim, name);
  struct exchanger* ex = mesh_ask_exchanger(m, dim);
  push_tag(ex, t, mesh_tags(m, dim));
}

void mesh_accumulate_tag(struct mesh* m, unsigned dim, const char* name)
{
  if (!mesh_is_parallel(m))
    return;
  struct const_tag* t = mesh_find_tag(m, dim, name);
  assert(t->type == TAG_F64);
  struct exchanger* ex = mesh_ask_exchanger(m, dim);
  double* out = exchange_doubles_add(ex, t->ncomps, t->d.f64,
      EX_REV, EX_ITEM);
  modify_tag(mesh_tags(m, dim), t->name, out);
}

unsigned mesh_ghost_layers(struct mesh* m)
{
  if (!mesh_is_parallel(m))
    return 0;
  return mesh_parallel(m)->nghost_layers;
}

void mesh_set_ghost_layers(struct mesh* m, unsigned n)
{
  mesh_parallel(m)->nghost_layers = n;
}

void mesh_set_globals(struct mesh* m, unsigned dim, unsigned long* new_globals)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  invalidate_globals(pm, dim);
  pm->globals[dim] = new_globals;
}

void mesh_set_own_ranks(struct mesh* m, unsigned dim, unsigned* new_owners)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  invalidate_ranks(pm, dim);
  pm->own_ranks[dim] = new_owners;
}

void mesh_tag_globals(struct mesh* m, unsigned dim)
{
  mesh_add_tag(m, dim, TAG_U64, "global_number", 1,
      copy_array(mesh_ask_globals(m, dim), mesh_count(m, dim)));
}

void mesh_tag_own_rank(struct mesh* m, unsigned dim)
{
  mesh_add_tag(m, dim, TAG_U32, "own_rank", 1,
      copy_array(mesh_ask_own_ranks(m, dim), mesh_count(m, dim)));
}

void mesh_tag_own_id(struct mesh* m, unsigned dim)
{
  mesh_add_tag(m, dim, TAG_U32, "own_id", 1,
      copy_array(mesh_ask_own_ids(m, dim), mesh_count(m, dim)));
}

void mesh_parallel_to_tags(struct mesh* m, unsigned dim)
{
  struct parallel_mesh* pm = mesh_parallel(m);
  if (pm->globals[dim])
    mesh_tag_globals(m, dim);
  if (pm->own_ranks[dim])
    mesh_tag_own_rank(m, dim);
  if (pm->own_ids[dim])
    mesh_tag_own_id(m, dim);
}

void mesh_parallel_untag(struct mesh* m, unsigned dim)
{
  if (mesh_find_tag(m, dim, "global_number"))
    mesh_free_tag(m, dim, "global_number");
  if (mesh_find_tag(m, dim, "own_rank"))
    mesh_free_tag(m, dim, "own_rank");
  if (mesh_find_tag(m, dim, "own_id"))
    mesh_free_tag(m, dim, "own_id");
}

void mesh_parallel_from_tags(struct mesh* m, unsigned dim)
{
  if (!mesh_find_tag(m, dim, "global_number"))
    return;
  mesh_set_globals(m, dim, copy_array(
        mesh_find_tag(m, dim, "global_number")->d.u64,
        mesh_count(m, dim)));
  mesh_free_tag(m, dim, "global_number");
  if (!mesh_find_tag(m, dim, "own_rank"))
    return;
  mesh_set_own_ranks(m, dim, copy_array(
        mesh_find_tag(m, dim, "own_rank")->d.u32,
        mesh_count(m, dim)));
  mesh_free_tag(m, dim, "own_rank");
  if (!mesh_find_tag(m, dim, "own_id"))
    return;
  struct parallel_mesh* pm = mesh_parallel(m);
  invalidate_ids(pm, dim);
  pm->own_ids[dim] = copy_array(
       mesh_find_tag(m, dim, "own_id")->d.u32,
       mesh_count(m, dim));
  mesh_free_tag(m, dim, "own_id");
}

void mesh_partition_out(struct mesh** p_m, unsigned is_source)
{
  *p_m = bcast_mesh_metadata(*p_m, is_source);
  if (is_source) {
    for (unsigned d = 0; d <= mesh_dim(*p_m); ++d)
      invalidate_ranks(mesh_parallel(*p_m), d);
  }
}

struct mesh* read_and_partition_serial_mesh(char const* filename)
{
  struct mesh* m = 0;
  if (comm_rank() == 0) {
    comm_use(comm_self());
    m = read_mesh_vtk(filename);
    assert(!mesh_is_parallel(m));
    mesh_make_parallel(m);
    comm_use(comm_world());
  }
  mesh_partition_out(&m, comm_rank() == 0);
  balance_mesh_inertial(m);
  return m;
}

LOOP_KERNEL(mark_owned,
    unsigned const* own_ranks,
    unsigned rank,
    unsigned* owned)
  owned[i] = (own_ranks[i] == rank);
}

unsigned* mesh_get_owned(struct mesh* m, unsigned dim)
{
  unsigned n = mesh_count(m, dim);
  unsigned* out = LOOP_MALLOC(unsigned, n);
  unsigned const* own_ranks = mesh_ask_own_ranks(m, dim);
  unsigned self = comm_rank();
  LOOP_EXEC(mark_owned, n, own_ranks, self, out);
  return out;
}

void mesh_doubles_max(struct mesh* m, unsigned dim, unsigned width,
    double** a)
{
  if (!mesh_is_parallel(m))
    return;
  double* in = *a;
  double* out = exchange_doubles_max(mesh_ask_exchanger(m, dim), width, in,
      EX_REV, EX_ITEM);
  loop_free(in);
  *a = out;
}

}
