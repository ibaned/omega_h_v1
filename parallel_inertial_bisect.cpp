#include "parallel_inertial_bisect.hpp"

#include <cassert>

#include "arrays.hpp"
#include "comm.hpp"
#include "doubles.hpp"
#include "element_field.hpp"
#include "exchanger.hpp"
#include "global.hpp"
#include "ghost_mesh.hpp"
#include "inertia.hpp"
#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "migrate_mesh.hpp"

void parallel_inertial_bisect(
    unsigned* p_n,
    double** p_coords,
    double** p_masses,
    unsigned** p_orig_ranks,
    unsigned** p_orig_ids)
{
  assert(comm_size() % 2 == 0);
  unsigned n = *p_n;
  double const* coords = *p_coords;
  double const* masses = 0;
  if (p_masses)
    masses = *p_masses;
  unsigned const* orig_ranks = *p_orig_ranks;
  unsigned const* orig_ids = *p_orig_ids;
  unsigned* marked = mark_inertial_bisection(n, coords, masses, 1);
  unsigned* offsets = uints_exscan(marked, n);
  loop_free(marked);
  unsigned nsubranks = comm_size() / 2;
  unsigned first_rank = 0;
  unsigned n_out = 0;
  double* coords_out = 0;
  double* masses_out = 0;
  unsigned* orig_ranks_out = 0;
  unsigned* orig_ids_out = 0;
  for (unsigned dir = 0; dir < 2; ++dir) {
    unsigned nsub = uints_at(offsets, n);
    unsigned* local = uints_linear(nsub + 1, 1);
    unsigned long* global = globalize_offsets(local, nsub);
    loop_free(local);
    unsigned* lin_ranks;
    unsigned* lin_ids;
    unsigned long nsub_total = comm_add_ulong(nsub);
    global_to_linpart(global, nsub, nsub_total, nsubranks,
        &lin_ranks, &lin_ids);
    loop_free(lin_ids);
    loop_free(global);
    for (unsigned i = 0; i < nsub; ++i)
      lin_ranks[i] += first_rank;
    unsigned in_subgroup = ((first_rank <= comm_rank()) &&
        (comm_rank() < (first_rank + nsubranks)));
    struct exchanger* ex = new_exchanger(nsub, lin_ranks);
    loop_free(lin_ranks);
    if (in_subgroup)
      n_out = ex->nitems[EX_REV];
    double* sub_coords = doubles_expand(n, 3, coords, offsets);
    double* coords_recvd = exchange_doubles(ex, 3, sub_coords, EX_FOR, EX_ITEM);
    loop_free(sub_coords);
    double* masses_recvd = 0;
    if (masses) {
      double* sub_masses = doubles_expand(n, 1, masses, offsets);
      masses_recvd = exchange_doubles(ex, 1, sub_masses, EX_FOR, EX_ITEM);
      loop_free(sub_masses);
    }
    unsigned* sub_orig_ranks = uints_expand(n, 1, orig_ranks, offsets);
    unsigned* orig_ranks_recvd = exchange_uints(ex, 1, sub_orig_ranks,
        EX_FOR, EX_ITEM);
    loop_free(sub_orig_ranks);
    unsigned* sub_orig_ids = uints_expand(n, 1, orig_ids, offsets);
    unsigned* orig_ids_recvd = exchange_uints(ex, 1, sub_orig_ids,
        EX_FOR, EX_ITEM);
    loop_free(sub_orig_ids);
    free_exchanger(ex);
    if (in_subgroup) {
      coords_out = coords_recvd;
      masses_out = masses_recvd;
      orig_ranks_out = orig_ranks_recvd;
      orig_ids_out = orig_ids_recvd;
    } else {
      loop_free(coords_recvd);
      loop_free(masses_recvd);
      loop_free(orig_ranks_recvd);
      loop_free(orig_ids_recvd);
    }
    unsigned* opp_offsets = uints_negate_offsets(offsets, n);
    loop_free(offsets);
    offsets = opp_offsets;
    first_rank += nsubranks;
  }
  loop_free(offsets);
  *p_n = n_out;
  loop_free(*p_coords);
  *p_coords = coords_out;
  if (p_masses) {
    loop_free(*p_masses);
    *p_masses = masses_out;
  }
  loop_free(*p_orig_ranks);
  *p_orig_ranks = orig_ranks_out;
  loop_free(*p_orig_ids);
  *p_orig_ids = orig_ids_out;
}

void recursive_inertial_bisect(
    unsigned* p_n,
    double** p_coords,
    double** p_masses,
    unsigned** p_orig_ranks,
    unsigned** p_orig_ids)
{
  if (comm_size() == 1)
    return;
  assert(comm_size() % 2 == 0);
  parallel_inertial_bisect(p_n, p_coords, p_masses,
      p_orig_ranks, p_orig_ids);
  unsigned half = comm_size() / 2;
  unsigned group = comm_rank() / half;
  unsigned subrank = comm_rank() % half;
  struct comm* oldcomm = comm_using();
  struct comm* subcomm = comm_split(oldcomm, group, subrank);
  comm_use(subcomm);
  recursive_inertial_bisect(p_n, p_coords, p_masses,
      p_orig_ranks, p_orig_ids);
  comm_use(oldcomm);
  comm_free(subcomm);
}

void balance_mesh_inertial(struct mesh* m)
{
  assert(mesh_is_parallel(m));
  mesh_ensure_ghosting(m, 0);
  unsigned dim = mesh_dim(m);
  unsigned had_elem_coords =
    (0 != mesh_find_tag(m, dim, "coordinates"));
  if (!had_elem_coords)
    mesh_interp_to_elems(m, "coordinates");
  unsigned n = mesh_count(m, dim);
  double* coords = generic_copy(
      mesh_find_tag(m, dim, "coordinates")->d.f64, n * 3);
  if (!had_elem_coords)
    mesh_free_tag(m, dim, "coordinates");
  unsigned* orig_ranks = uints_filled(n, comm_rank());
  unsigned* orig_ids = uints_linear(n, 1);
  recursive_inertial_bisect(&n, &coords, 0, &orig_ranks, &orig_ids);
  loop_free(coords);
  migrate_mesh(m, n, orig_ranks, orig_ids);
  loop_free(orig_ranks);
  loop_free(orig_ids);
}
