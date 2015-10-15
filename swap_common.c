#include "swap_common.h"

#include <assert.h>

#include "concat.h"
#include "doubles.h"
#include "graph.h"
#include "indset.h"
#include "ints.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "quality.h"
#include "swap_qualities.h"
#include "swap_topology.h"
#include "tag.h"

unsigned swap_common(
    struct mesh** p_m,
    unsigned* candidates)
{
  struct mesh* m = *p_m;
  unsigned elem_dim = mesh_dim(m);
  assert(elem_dim == 3);
  unsigned nedges = mesh_count(m, 1);
  if (!uints_max(candidates, nedges))
    return 0;
  unsigned const* verts_of_edges = mesh_ask_down(m, 1, 0);
  unsigned const* class_dim = mesh_find_tag(m, 0, "class_dim")->d.u32;
  unmark_boundary(elem_dim, 1, nedges, verts_of_edges, class_dim, candidates);
  if (!uints_max(candidates, nedges))
    return 0;
  unsigned const* elems_of_edges_offsets =
    mesh_ask_up(m, 1, elem_dim)->offsets;
  unsigned const* elems_of_edges =
    mesh_ask_up(m, 1, elem_dim)->adj;
  unsigned const* elems_of_edges_directions =
    mesh_ask_up(m, 1, elem_dim)->directions;
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  double const* coords = mesh_find_tag(m, 0, "coordinates")->d.f64;
  double* elem_quals = mesh_qualities(m);
  double* edge_quals;
  unsigned* edge_codes;
  unsigned* gen_elems_per_edge;
  swap_qualities(nedges, candidates,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions,
      verts_of_edges, verts_of_elems, coords, elem_quals,
      &edge_quals, &edge_codes, &gen_elems_per_edge);
  loop_free(elem_quals);
  if (!uints_max(candidates, nedges)) {
    loop_free(edge_quals);
    loop_free(edge_codes);
    loop_free(gen_elems_per_edge);
    return 0;
  }
  unsigned const* edges_of_edges_offsets = mesh_ask_star(m, 1, elem_dim)->offsets;
  unsigned const* edges_of_edges = mesh_ask_star(m, 1, elem_dim)->adj;
  unsigned* indset = find_indset(nedges, edges_of_edges_offsets, edges_of_edges,
      candidates, edge_quals);
  loop_free(edge_quals);
  for (unsigned i = 0; i < nedges; ++i)
    if (!indset[i])
      gen_elems_per_edge[i] = 0;
  unsigned* gen_offset_of_edges = uints_exscan(gen_elems_per_edge, nedges);
  unsigned ngen_elems = gen_offset_of_edges[nedges];
  loop_free(gen_elems_per_edge);
  unsigned* verts_of_gen_elems = swap_topology(nedges, indset,
      gen_offset_of_edges, edge_codes,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions,
      verts_of_edges, verts_of_elems);
  loop_free(edge_codes);
  loop_free(gen_offset_of_edges);
  unsigned* old_elems = mesh_mark_up(m, 1, elem_dim, indset);
  loop_free(indset);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned* same_elems = uints_negate(old_elems, nelems);
  loop_free(old_elems);
  unsigned* same_elem_offsets = uints_exscan(same_elems, nelems);
  loop_free(same_elems);
  unsigned nelems_out;
  unsigned* verts_of_elems_out;
  concat_verts_of_elems(elem_dim, nelems, ngen_elems, verts_of_elems,
      same_elem_offsets, verts_of_gen_elems,
      &nelems_out, &verts_of_elems_out);
  loop_free(same_elem_offsets);
  loop_free(verts_of_gen_elems);
  unsigned nverts = mesh_count(m, 0);
  struct mesh* m_out = new_mesh(elem_dim);
  mesh_set_ents(m_out, 0, nverts, 0);
  mesh_set_ents(m_out, elem_dim, nelems_out, verts_of_elems_out);
  for (unsigned i = 0; i < mesh_count_tags(m, 0); ++i) {
    struct const_tag* t = mesh_get_tag(m, 0, i);
    void* data = 0;
    switch (t->type) {
      case TAG_U8:  data = uchars_copy(t->d.u8, nverts * t->ncomps);
                    break;
      case TAG_U32: data = uints_copy(t->d.u32, nverts * t->ncomps);
                    break;
      case TAG_U64: data = ulongs_copy(t->d.u64, nverts * t->ncomps);
                    break;
      case TAG_F64: data = doubles_copy(t->d.f64, nverts * t->ncomps);
                    break;
    };
    mesh_add_tag(m_out, 0, t->type, t->name, t->ncomps, data);
  }
  free_mesh(m);
  *p_m = m_out;
  return 1;
}
