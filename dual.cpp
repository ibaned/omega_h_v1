#include "dual.hpp"

#include <cassert>

#include "ints.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tables.hpp"

LOOP_INOUT static unsigned copy_except(
    unsigned const* a,
    unsigned* b,
    unsigned n,
    unsigned exclude_this)
{
  unsigned j = 0;
  for (unsigned i = 0; i < n; ++i)
    if (a[i] != exclude_this)
      b[j++] = a[i];
  return j;
}

LOOP_INOUT static unsigned intersect(
    unsigned* a,
    unsigned na,
    unsigned const* b,
    unsigned nb)
{
  unsigned j = 0;
  for (unsigned i = 0; i < na; ++i)
    if (has(b, nb, a[i]))
      a[j++] = a[i];
  return j;
}

LOOP_KERNEL(element_dual_from_verts,
    unsigned elem_dim,
    unsigned side_dim,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    unsigned verts_per_elem,
    unsigned sides_per_elem,
    unsigned verts_per_side,
    unsigned* elems_of_elems)

  unsigned const* const* elem_verts_of_sides =
      the_canonical_orders[elem_dim][side_dim][0];
  unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
  unsigned* elems_of_elem = elems_of_elems + i * sides_per_elem;
  for (unsigned j = 0; j < sides_per_elem; ++j) {
    unsigned const* elem_verts_of_side = elem_verts_of_sides[j];
    unsigned buf[MAX_UP];
    unsigned buf_size = 0;
    for (unsigned k = 0; k < verts_per_side; ++k) {
      unsigned vert = verts_of_elem[elem_verts_of_side[k]];
      unsigned first_use = elems_of_verts_offsets[vert];
      unsigned end_use = elems_of_verts_offsets[vert + 1];
      if (k) {
        buf_size = intersect(
            buf,
            buf_size,
            elems_of_verts + first_use,
            end_use - first_use);
      } else {
        assert(end_use - first_use <= MAX_UP);
        buf_size = copy_except(
            elems_of_verts + first_use,
            buf,
            end_use - first_use,
            i);
      }
    }
    assert(buf_size <= 1);
    elems_of_elem[j] = ((buf_size) ? buf[0] : INVALID);
  }
}

unsigned* dual_from_verts(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts)
{
  unsigned side_dim = elem_dim - 1;
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned sides_per_elem = the_down_degrees[elem_dim][side_dim];
  unsigned verts_per_side = the_down_degrees[side_dim][0];
  unsigned* elems_of_elems = LOOP_MALLOC(unsigned, nelems * sides_per_elem);
  LOOP_EXEC(element_dual_from_verts, nelems,
      elem_dim,
      side_dim,
      verts_of_elems,
      elems_of_verts_offsets,
      elems_of_verts,
      verts_per_elem,
      sides_per_elem,
      verts_per_side,
      elems_of_elems);
  return elems_of_elems;
}

unsigned* mesh_get_dual_from_verts(struct mesh* m)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  struct const_up* elems_of_verts = mesh_ask_up(m, 0, elem_dim);
  return dual_from_verts(elem_dim, nelems, verts_of_elems,
      elems_of_verts->offsets, elems_of_verts->adj);
}

LOOP_KERNEL(element_dual_from_sides,
    unsigned sides_per_elem,
    unsigned const* sides_of_elems,
    unsigned const* elems_of_sides_offsets,
    unsigned const* elems_of_sides,
    unsigned* elems_of_elems)

  unsigned const* sides_of_elem = sides_of_elems + i * sides_per_elem;
  unsigned* elems_of_elem = elems_of_elems + i * sides_per_elem;
  for (unsigned j = 0; j < sides_per_elem; ++j) {
    unsigned side = sides_of_elem[j];
    unsigned f = elems_of_sides_offsets[side];
    unsigned e = elems_of_sides_offsets[side + 1];
    unsigned k;
    for (k = f; k < e; ++k) {
      unsigned oelem = elems_of_sides[k];
      if (oelem != i) {
        elems_of_elem[j] = oelem;
        break;
      }
    }
    if (k == e)
      elems_of_elem[j] = INVALID;
  }
}

unsigned* dual_from_sides(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* sides_of_elems,
    unsigned const* elems_of_sides_offsets,
    unsigned const* elems_of_sides)
{
  unsigned sides_per_elem = the_down_degrees[elem_dim][elem_dim - 1];
  unsigned* elems_of_elems = LOOP_MALLOC(unsigned, nelems * sides_per_elem);
  LOOP_EXEC(element_dual_from_sides, nelems,
      sides_per_elem,
      sides_of_elems,
      elems_of_sides_offsets,
      elems_of_sides,
      elems_of_elems);
  return elems_of_elems;
}

unsigned* mesh_get_dual_from_sides(struct mesh* m)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* sides_of_elems = mesh_ask_down(m, elem_dim, elem_dim - 1);
  struct const_up* elems_of_sides = mesh_ask_up(m, elem_dim - 1, elem_dim);
  return dual_from_sides(elem_dim, nelems, sides_of_elems,
      elems_of_sides->offsets, elems_of_sides->adj);
}
