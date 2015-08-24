#include "refine_common.h"
#include "refine_qualities.h"
#include "indset.h"
#include "ints.h"
#include "splits_to_elements.h"
#include "refine_topology.h"
#include "refine_nodal.h"
#include "refine_class.h"
#include "concat.h"
#include <stdlib.h>

void refine_common(
    struct mesh** p_m,
    unsigned src_dim,
    unsigned const* candidates)
{
  struct mesh* m = *p_m;
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
  double const* coords = mesh_find_nodal_field(m, "coordinates")->data;
  double* src_quals = refine_qualities(elem_dim, src_dim, nsrcs, verts_of_srcs,
      verts_of_elems, elems_of_srcs_offsets, elems_of_srcs,
      elems_of_srcs_directions, candidates, coords);
  unsigned const* srcs_of_srcs_offsets =
    mesh_ask_star(m, src_dim, elem_dim)->offsets;
  unsigned const* srcs_of_srcs =
    mesh_ask_star(m, src_dim, elem_dim)->adj;
  unsigned* indset = find_indset(nsrcs, srcs_of_srcs_offsets, srcs_of_srcs,
      candidates, src_quals);
  free(src_quals);
  unsigned* gen_offset_of_srcs = ints_exscan(indset, nsrcs);
  free(indset);
  unsigned nsplit_srcs = gen_offset_of_srcs[nsrcs];
  unsigned nverts = mesh_count(m, 0);
  unsigned* gen_vert_of_srcs = malloc(sizeof(unsigned) * nsrcs);
  for (unsigned i = 0; i < nsrcs; ++i)
    if (gen_offset_of_srcs[i] != gen_offset_of_srcs[i + src_dim])
      gen_vert_of_srcs[i] = nverts + gen_offset_of_srcs[i];
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* srcs_of_elems = mesh_ask_down(m, elem_dim, src_dim);
  unsigned* gen_offset_of_elems;
  unsigned* gen_direction_of_elems;
  unsigned* gen_vert_of_elems;
  project_splits_to_elements(elem_dim, src_dim, nelems,
      srcs_of_elems, gen_offset_of_srcs, gen_vert_of_srcs,
      &gen_offset_of_elems, &gen_direction_of_elems, &gen_vert_of_elems);
  free(gen_vert_of_srcs);
  unsigned ngen_elems;
  unsigned* verts_of_gen_elems;
  refine_topology(elem_dim, src_dim, elem_dim, nelems, verts_of_elems,
      gen_offset_of_elems, gen_vert_of_elems, gen_direction_of_elems,
      &ngen_elems, &verts_of_gen_elems);
  free(gen_vert_of_elems);
  free(gen_direction_of_elems);
  struct mesh* m_out = new_mesh(elem_dim);
  unsigned nverts_out = nverts + nsplit_srcs;
  mesh_set_ents(m_out, 0, nverts_out, 0);
  for (unsigned i = 0; i < mesh_count_nodal_fields(m); ++i) {
    struct const_field* f = mesh_get_nodal_field(m, i);
    double* gen_vals = refine_nodal(src_dim, nsrcs, verts_of_srcs,
        gen_offset_of_srcs, f->ncomps, f->data);
    double* vals_out = concat_doubles(f->ncomps, f->data, nverts,
        gen_vals, nsplit_srcs);
    free(gen_vals);
    mesh_add_nodal_field(m_out, f->name, f->ncomps, vals_out);
  }
  if (mesh_find_nodal_label(m, "class_dim")) {
    unsigned const* class_dim = mesh_find_nodal_label(m, "class_dim")->data;
    unsigned* gen_class_dim = refine_class(src_dim, nsrcs, verts_of_srcs,
        gen_offset_of_srcs, class_dim);
    unsigned* class_dim_out = concat_ints(1, class_dim, nverts,
        gen_class_dim, nsplit_srcs);
    free(gen_class_dim);
    mesh_add_nodal_label(m_out, "class_dim", class_dim_out);
  }
  free(gen_offset_of_srcs);
  unsigned* offset_of_same_elems = ints_negate_offsets(
      gen_offset_of_elems, nelems);
  free(gen_offset_of_elems);
  unsigned nelems_out;
  unsigned* verts_of_elems_out;
  concat_verts_of_elems(elem_dim, nelems, ngen_elems, verts_of_elems,
      offset_of_same_elems, verts_of_gen_elems,
      &nelems_out, &verts_of_elems_out);
  mesh_set_ents(m_out, elem_dim, nelems_out, verts_of_elems_out);
  free(offset_of_same_elems);
  free(verts_of_gen_elems);
  free_mesh(m);
  *p_m = m_out;
}
