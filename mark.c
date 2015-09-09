#include "mark.h"
#include "ints.h"
#include "mesh.h"
#include "tables.h"
#include <stdlib.h>

unsigned* mark_down(
    unsigned nlows,
    unsigned const* highs_of_lows_offsets,
    unsigned const* highs_of_lows,
    unsigned const* marked_highs)
{
  unsigned* low_marks = malloc(sizeof(unsigned) * nlows);
  for (unsigned i = 0; i < nlows; ++i) {
    unsigned first_use = highs_of_lows_offsets[i];
    unsigned end_use = highs_of_lows_offsets[i + 1];
    low_marks[i] = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned high = highs_of_lows[j];
      if (marked_highs[high]) {
        low_marks[i] = 1;
        break;
      }
    }
  }
  return low_marks;
}

unsigned* mesh_mark_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned const* marked_highs)
{
  return mark_down(mesh_count(m, low_dim),
      mesh_ask_up(m, low_dim, high_dim)->offsets,
      mesh_ask_up(m, low_dim, high_dim)->adj,
      marked_highs);
}

static unsigned* mark_dual(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* dual,
    unsigned const* marked)
{
  unsigned degree = the_down_degrees[elem_dim][elem_dim - 1];
  unsigned* out = malloc(sizeof(unsigned) * nelems);
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
    free(in);
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

void unmark_boundary(
    unsigned elem_dim,
    unsigned ent_dim,
    unsigned nents,
    unsigned const* verts_of_ents,
    unsigned const* vert_class_dim,
    unsigned* marked)
{
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  for (unsigned i = 0; i < nents; ++i) {
    if (!marked[i])
      continue;
    unsigned const* verts_of_ent = verts_of_ents + i * verts_per_ent;
    unsigned is_boundary = 1;
    for (unsigned j = 0; j < verts_per_ent; ++j) {
      unsigned vert = verts_of_ent[j];
      if (vert_class_dim[vert] == elem_dim)
        is_boundary = 0;
    }
    if (is_boundary)
      marked[i] = 0;
  }
}
