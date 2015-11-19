#include "mark.h"

#include <assert.h>

#include "infer_class.h"
#include "loop.h"
#include "mesh.h"
#include "quality.h"
#include "tables.h"

unsigned* mark_down(
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* marked_highs)
{
  unsigned* marked_lows = LOOP_MALLOC(unsigned, nlows);
  for (unsigned i = 0; i < nlows; ++i) {
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

unsigned* mesh_mark_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* marked_highs)
{
  return mark_down(mesh_count(m, low_dim),
      mesh_ask_up(m, low_dim, high_dim)->offsets,
      mesh_ask_up(m, low_dim, high_dim)->adj,
      marked_highs);
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

static void mark_dual_layers(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* dual,
    unsigned** marked,
    unsigned nlayers)
{
  for (unsigned i = 0; i < nlayers; ++i) {
    unsigned* in = *marked;
    unsigned* out = mark_dual(elem_dim, nelems, dual, in);
    loop_free(in);
    *marked = out;
  }
}

void mesh_mark_dual_layers(
    struct mesh* m,
    unsigned** marked,
    unsigned nlayers)
{
  mark_dual_layers(mesh_dim(m), mesh_count(m, mesh_dim(m)),
      mesh_ask_dual(m), marked, nlayers);
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
  for (unsigned i = 0; i < nents; ++i) {
    out[i] = 0;
    if (class_dim_of_ents[i] == target_dim)
      if (target_id == INVALID || class_id_of_ents[i] == target_id)
        out[i] = 1;
  }
  return out;
}

unsigned* mesh_mark_class(struct mesh* m, unsigned ent_dim,
    unsigned target_dim, unsigned target_id)
{
  unsigned const* class_dim_of_ents = mesh_ask_class_dim(m, ent_dim);
  unsigned const* class_id_of_ents = mesh_ask_class_id(m, ent_dim);
  unsigned nents = mesh_count(m, ent_dim);
  return mark_class(nents, target_dim, target_id, class_dim_of_ents,
      class_id_of_ents);
}

unsigned* mesh_mark_class_closure_verts(struct mesh* m, unsigned target_dim,
    unsigned target_id)
{
  unsigned* equal_order = mesh_mark_class(m, target_dim, target_dim, target_id);
  return mesh_mark_down(m, target_dim, 0, equal_order);
}

void unmark_boundary(
    unsigned elem_dim,
    unsigned ent_dim,
    unsigned nents,
    unsigned const* verts_of_ents,
    unsigned const* vert_class_dim,
    unsigned* marked)
{
  unsigned* ent_class_dim;
  infer_class(ent_dim, nents, verts_of_ents, vert_class_dim, 0,
      &ent_class_dim, 0);
  unsigned* interior_ents = mark_class(nents, elem_dim, INVALID,
      ent_class_dim, 0);
  for (unsigned i = 0; i < nents; ++i)
    if (!interior_ents[i])
      marked[i] = 0;
}

static unsigned* mark_slivers(
    unsigned nelems,
    double const* elem_quals,
    double good_qual)
{
  unsigned* slivers = LOOP_MALLOC(unsigned, nelems);
  for (unsigned i = 0; i < nelems; ++i)
    slivers[i] = (elem_quals[i] < good_qual) ? 1 : 0;
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

unsigned* mark_part_boundary(
    unsigned nsides,
    unsigned const* elems_of_sides_offsets)
{
  unsigned* out = LOOP_MALLOC(unsigned, nsides);
  for (unsigned i = 0; i < nsides; ++i)
    out[i] = ( elems_of_sides_offsets[i + 1]
             - elems_of_sides_offsets[i]    ) < 2;
  return out;
}

unsigned* mesh_mark_part_boundary(struct mesh* m)
{
  unsigned dim = mesh_dim(m);
  return mark_part_boundary(mesh_count(m, dim - 1),
      mesh_ask_up(m, dim - 1, dim)->offsets);
}
