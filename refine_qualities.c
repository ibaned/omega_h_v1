#include "refine_qualities.h"

#include <assert.h>

#include "algebra.h"
#include "element_field.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "quality.h"
#include "tables.h"

double* refine_qualities(
    unsigned elem_dim,
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_srcs_offsets,
    unsigned const* elems_of_srcs,
    unsigned const* elems_of_srcs_directions,
    unsigned* candidate_srcs,
    double const* coords,
    double qual_floor,
    double const* elem_quals,
    unsigned require_better)
{
  assert(elem_dim >= src_dim);
  assert(src_dim > 0);
  unsigned base_dim = elem_dim - 1;
  unsigned opp_dim = get_opposite_dim(elem_dim, base_dim);
  assert(opp_dim == 0);
  unsigned verts_per_src = the_down_degrees[src_dim][0];
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  double* out = LOOP_MALLOC(double, nsrcs);
  unsigned const* const* elem_verts_of_srcs =
    the_canonical_orders[elem_dim][src_dim][0];
  unsigned const* const* elem_verts_of_bases =
    the_canonical_orders[elem_dim][base_dim][0];
  unsigned const* elem_base_of_opps = the_opposite_orders[elem_dim][0];
  quality_function qf = the_quality_functions[elem_dim];
  for (unsigned i = 0; i < nsrcs; ++i) {
    if (!candidate_srcs[i])
      continue;
    unsigned first_use = elems_of_srcs_offsets[i];
    unsigned end_use = elems_of_srcs_offsets[i + 1];
    unsigned const* verts_of_src = verts_of_srcs + i * verts_per_src;
    double split_x[3];
    average_element_field(verts_per_src, verts_of_src,
        coords, 3, split_x);
    double minq = 1;
    double old_minq = 1;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned direction = elems_of_srcs_directions[j];
      unsigned elem = elems_of_srcs[j];
      if (require_better) {
        double old_q = elem_quals[elem];
        if (old_q < old_minq)
          old_minq = old_q;
      }
      unsigned const* verts_of_elem = verts_of_elems + elem * verts_per_elem;
      unsigned const* elem_verts_of_src = elem_verts_of_srcs[direction];
      for (unsigned k = 0; k < verts_per_src; ++k) {
        unsigned opp = elem_verts_of_src[k];
        unsigned base = elem_base_of_opps[opp];
        unsigned const* elem_verts_of_base = elem_verts_of_bases[base];
        double elem_x[MAX_DOWN][3];
        for (unsigned l = 0; l < (verts_per_elem - 1); ++l) {
          unsigned vert = verts_of_elem[elem_verts_of_base[l]];
          copy_vector(coords + vert * 3, elem_x[l], 3);
        }
        copy_vector(split_x, elem_x[verts_per_elem - 1], 3);
        double q = qf(elem_x);
        assert(q > 0);
        if (q < minq)
          minq = q;
      }
    }
    if ((minq < qual_floor) ||
        (require_better && (minq <= old_minq)))
      candidate_srcs[i] = 0;
    else
      out[i] = minq;
  }
  return out;
}

double* mesh_refine_qualities(struct mesh* m, unsigned src_dim,
    unsigned* candidates, double qual_floor, unsigned require_better)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nsrcs = mesh_count(m, src_dim);
  unsigned const* verts_of_srcs = mesh_ask_down(m, src_dim, 0);
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned const* elems_of_srcs_offsets =
    mesh_ask_up(m, src_dim, elem_dim)->offsets;
  unsigned const* elems_of_srcs =
    mesh_ask_up(m, src_dim, elem_dim)->adj;
  unsigned const* elems_of_srcs_directions =
    mesh_ask_up(m, src_dim, elem_dim)->directions;
  double const* coords = mesh_find_tag(m, 0, "coordinates")->d.f64;
  double* elem_quals = 0;
  if (require_better)
    elem_quals = mesh_qualities(m);
  double* src_quals = refine_qualities(elem_dim, src_dim, nsrcs,
      verts_of_srcs, verts_of_elems,
      elems_of_srcs_offsets, elems_of_srcs, elems_of_srcs_directions,
      candidates, coords, qual_floor,
      elem_quals, require_better);
  loop_free(elem_quals);
  if (mesh_is_parallel(m))
    mesh_conform_doubles(m, src_dim, 1, &src_quals);
  return src_quals;
}
