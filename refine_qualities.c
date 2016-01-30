#include "refine_qualities.h"

#include <assert.h>

#include "algebra.h"
#include "element_field.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "quality.h"
#include "tables.h"

LOOP_KERNEL(refine_quality,
    unsigned const* elems_of_srcs_offsets,
    unsigned const* elems_of_srcs,
    unsigned const* elems_of_srcs_directions,
    unsigned const* verts_of_srcs,
    unsigned verts_per_src,
    double const* coords,
    double const* elem_quals,
    unsigned const* verts_of_elems,
    unsigned verts_per_elem,
    unsigned** elem_verts_of_srcs,
    unsigned** elem_verts_of_bases,
    unsigned* elem_base_of_opps,
    quality_function qf,
    double qual_floor,
    unsigned* candidate_srcs,
    double* out)

  if (!candidate_srcs[i])
    return;
  unsigned first_use = elems_of_srcs_offsets[i];
  unsigned end_use = elems_of_srcs_offsets[i + 1];
  unsigned const* verts_of_src = verts_of_srcs + i * verts_per_src;
  double split_x[3];
  average_element_field(verts_per_src, verts_of_src,
      coords, 3, split_x);
  double minq = 1;
  double old_minq = 1;
  unsigned require_better = (elem_quals != 0);
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

double* mesh_refine_qualities(struct mesh* m, unsigned src_dim,
    unsigned** p_candidates, double qual_floor, unsigned require_better)
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
  assert(elem_dim >= src_dim);
  assert(src_dim > 0);
  unsigned base_dim = elem_dim - 1;
  unsigned opp_dim = get_opposite_dim(elem_dim, base_dim);
  assert(opp_dim == 0);
  unsigned verts_per_src = the_down_degrees[src_dim][0];
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  double* src_quals = LOOP_MALLOC(double, nsrcs);
  unsigned** elem_verts_of_srcs = orders_to_device(elem_dim, src_dim, 0);
  unsigned** elem_verts_of_bases = orders_to_device(elem_dim, base_dim, 0);
  unsigned* elem_base_of_opps = LOOP_TO_DEVICE(unsigned,
      the_opposite_orders[elem_dim][opp_dim], verts_per_elem);
  quality_function qf = the_quality_functions[elem_dim];
  LOOP_EXEC(refine_quality, nsrcs,
      elems_of_srcs_offsets,
      elems_of_srcs,
      elems_of_srcs_directions,
      verts_of_srcs,
      verts_per_src,
      coords,
      elem_quals,
      verts_of_elems,
      verts_per_elem,
      elem_verts_of_srcs,
      elem_verts_of_bases,
      elem_base_of_opps,
      qf,
      qual_floor,
      *p_candidates,
      src_quals);
  free_orders(elem_verts_of_srcs, elem_dim, src_dim);
  free_orders(elem_verts_of_bases, elem_dim, base_dim);
  loop_host_free(elem_base_of_opps);
  loop_free(elem_quals);
  mesh_conform_doubles(m, src_dim, 1, &src_quals);
  mesh_conform_uints(m, src_dim, 1, p_candidates);
  return src_quals;
}
