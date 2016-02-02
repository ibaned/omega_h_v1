#include "mark.h"

#include <assert.h>

#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "quality.h"
#include "tables.h"

LOOP_KERNEL(mark_down_kern,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* marked_highs,
    unsigned* marked_lows)
  unsigned first_use = highs_of_lows_offsets[i];
  unsigned end_use = highs_of_lows_offsets[i + 1];
  marked_lows[i] = 0;
  for (unsigned j = first_use; j < end_use; ++j) {
    unsigned high = highs_of_lows[j];
    if (marked_highs[high]) {
      marked_lows[i] = 1;
      break;
    }
  }
}

unsigned* mark_down(
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* marked_highs)
{
  unsigned* marked_lows = LOOP_MALLOC(unsigned, nlows);
  LOOP_EXEC(mark_down_kern, nlows,
      highs_of_lows_offsets, highs_of_lows,
      marked_highs, marked_lows);
  return marked_lows;
}

unsigned* mark_up(
    unsigned high_dim,
    unsigned low_dim,
    unsigned nhighs,
    unsigned const* lows_of_highs,
    unsigned const* marked_lows)
{
  unsigned lows_per_high = the_down_degrees[high_dim][low_dim];
  unsigned* marked_highs = LOOP_MALLOC(unsigned, nhighs);
  for (unsigned i = 0; i < nhighs; ++i) {
    marked_highs[i] = 0;
    unsigned const* lows_of_high = lows_of_highs + i * lows_per_high;
    for (unsigned j = 0; j < lows_per_high; ++j)
      if (marked_lows[lows_of_high[j]]) {
        marked_highs[i] = 1;
        break;
      }
  }
  return marked_highs;
}

unsigned* mesh_mark_down_local(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* marked_highs)
{
  unsigned* out = mark_down(mesh_count(m, low_dim),
      mesh_ask_up(m, low_dim, high_dim)->offsets,
      mesh_ask_up(m, low_dim, high_dim)->adj,
      marked_highs);
  return out;
}

unsigned* mesh_mark_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* marked_highs)
{
  unsigned* out = mesh_mark_down_local(m, high_dim, low_dim, marked_highs);
  if (mesh_is_parallel(m))
    assert(mesh_ghost_layers(m) == 1);
  mesh_conform_uints(m, low_dim, 1, &out);
  return out;
}

unsigned* mesh_mark_up(struct mesh* m, unsigned low_dim, unsigned high_dim,
    unsigned const* marked_lows)
{
  return mark_up(high_dim, low_dim, mesh_count(m, high_dim),
      mesh_ask_down(m, high_dim, low_dim), marked_lows);
}

static unsigned* mark_dual(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* dual,
    unsigned const* marked)
{
  unsigned degree = the_down_degrees[elem_dim][elem_dim - 1];
  unsigned* out = LOOP_MALLOC(unsigned, nelems);
  for (unsigned i = 0; i < nelems; ++i) {
    out[i] = marked[i];
    unsigned const* elems_of_elem = dual + i * degree;
    for (unsigned j = 0; j < degree; ++j)
      if (elems_of_elem[j] != INVALID && marked[elems_of_elem[j]])
        out[i] = 1;
  }
  return out;
}

void mesh_mark_dual_layers(
    struct mesh* m,
    unsigned** marked,
    unsigned nlayers)
{
  if (mesh_is_parallel(m))
    assert(mesh_ghost_layers(m) == 1);
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* dual = mesh_ask_dual(m);
  for (unsigned i = 0; i < nlayers; ++i) {
    unsigned* in = *marked;
    unsigned* out = mark_dual(elem_dim, nelems, dual, in);
    loop_free(in);
    if (mesh_is_parallel(m))
      mesh_conform_uints(m, elem_dim, 1, &out);
    *marked = out;
  }
}

LOOP_KERNEL(mark_class_kern,
    unsigned target_dim,
    unsigned target_id,
    unsigned const* class_dim_of_ents,
    unsigned const* class_id_of_ents,
    unsigned* out)
  out[i] = 0;
  if (class_dim_of_ents[i] == target_dim)
    if (target_id == INVALID || class_id_of_ents[i] == target_id)
      out[i] = 1;
}

unsigned* mark_class(
    unsigned nents,
    unsigned target_dim,
    unsigned target_id,
    unsigned const* class_dim_of_ents,
    unsigned const* class_id_of_ents)
{
  if (target_id != INVALID)
    assert(class_id_of_ents);
  unsigned* out = LOOP_MALLOC(unsigned, nents);
  LOOP_EXEC(mark_class_kern, nents,
      target_dim, target_id, class_dim_of_ents, class_id_of_ents, out);
  return out;
}

unsigned* mesh_mark_class(struct mesh* m, unsigned ent_dim,
    unsigned target_dim, unsigned target_id)
{
  unsigned const* class_dim_of_ents =
    mesh_find_tag(m, ent_dim, "class_dim")->d.u32;
  unsigned const* class_id_of_ents =
    mesh_find_tag(m, ent_dim, "class_id")->d.u32;
  unsigned nents = mesh_count(m, ent_dim);
  return mark_class(nents, target_dim, target_id, class_dim_of_ents,
      class_id_of_ents);
}

unsigned* mesh_mark_class_closure_verts(struct mesh* m, unsigned target_dim,
    unsigned target_id)
{
  unsigned* equal_order = mesh_mark_class(m, target_dim, target_dim, target_id);
  unsigned* out = mesh_mark_down(m, target_dim, 0, equal_order);
  loop_free(equal_order);
  return out;
}

void mesh_unmark_boundary(
    struct mesh* m,
    unsigned ent_dim,
    unsigned* marked)
{
  unsigned const* class_dim =
    mesh_find_tag(m, ent_dim, "class_dim")->d.u32;
  unsigned elem_dim = mesh_dim(m);
  unsigned nents = mesh_count(m, ent_dim);
  for (unsigned i = 0; i < nents; ++i)
    if (class_dim[i] != elem_dim)
      marked[i] = 0;
}

LOOP_KERNEL(mark_sliver,
    double const* elem_quals,
    double good_qual,
    unsigned* slivers)
  slivers[i] = (elem_quals[i] < good_qual) ? 1 : 0;
}

static unsigned* mark_slivers(
    unsigned nelems,
    double const* elem_quals,
    double good_qual)
{
  unsigned* slivers = LOOP_MALLOC(unsigned, nelems);
  LOOP_EXEC(mark_sliver, nelems, elem_quals, good_qual, slivers);
  return slivers;
}

unsigned* mesh_mark_slivers(struct mesh* m, double good_qual, unsigned nlayers)
{
  unsigned nelems = mesh_count(m, mesh_dim(m));
  double* elem_quals = mesh_qualities(m);
  unsigned* slivers = mark_slivers(nelems, elem_quals, good_qual);
  loop_free(elem_quals);
  mesh_mark_dual_layers(m, &slivers, nlayers);
  return slivers;
}

LOOP_KERNEL(part_boundary_kern,
    unsigned const* elems_of_sides_offsets,
    unsigned* out)
  out[i] = ( elems_of_sides_offsets[i + 1]
           - elems_of_sides_offsets[i]    ) < 2;
}

unsigned* mark_part_boundary(
    unsigned nsides,
    unsigned const* elems_of_sides_offsets)
{
  unsigned* out = LOOP_MALLOC(unsigned, nsides);
  LOOP_EXEC(part_boundary_kern, nsides, elems_of_sides_offsets, out);
  return out;
}

unsigned* mesh_mark_part_boundary(struct mesh* m)
{
  unsigned dim = mesh_dim(m);
  return mark_part_boundary(mesh_count(m, dim - 1),
      mesh_ask_up(m, dim - 1, dim)->offsets);
}
