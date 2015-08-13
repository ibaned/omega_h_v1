#include "split_sliver_tris.h"
#include "bad_elem_keys.h"
#include "ints.h"
#include "up_from_down.h"
#include "star.h"
#include "bridge_graph.h"
#include "reflect_down.h"
#include "collect_keys.h"
#include "refine_common.h"
#include <assert.h>
#include <stdlib.h>

#include <stdio.h>

int split_sliver_tris(
    unsigned elem_dim,
    unsigned* p_nelems,
    unsigned* p_nverts,
    unsigned** p_verts_of_elems,
    double** p_coords,
    double qual_floor,
    double edge_ratio_floor)
{
  /* TODO: accept a tet mesh */
  assert(elem_dim == 2);
  unsigned nelems = *p_nelems;
  unsigned nverts = *p_nverts;
  unsigned const* verts_of_elems = *p_verts_of_elems;
  unsigned ntris = *p_nelems;
  unsigned const* verts_of_tris = *p_verts_of_elems;
  double const* coords = *p_coords;
  unsigned* bad_tris;
  unsigned* key_of_tris;
  bad_elem_keys(2, ntris, verts_of_tris, coords,
      SLIVER_ELEM, qual_floor, edge_ratio_floor,
      &bad_tris, &key_of_tris);
  unsigned something_to_do = ints_max(bad_tris, ntris);
  if (!something_to_do) {
    printf("no sliver tris found\n");
    free(bad_tris);
    free(key_of_tris);
    return 0;
  }
  unsigned* elems_of_verts_offsets;
  unsigned* elems_of_verts;
  up_from_down(elem_dim, 0, nelems, nverts, verts_of_elems,
      &elems_of_verts_offsets, &elems_of_verts, 0);
  unsigned* verts_of_verts_offsets;
  unsigned* verts_of_verts;
  get_star(0, elem_dim, nverts, elems_of_verts_offsets, elems_of_verts,
      verts_of_elems, &verts_of_verts_offsets, &verts_of_verts);
  free(elems_of_verts_offsets);
  free(elems_of_verts);
  unsigned nedges;
  unsigned* verts_of_edges;
  bridge_graph(nverts, verts_of_verts_offsets, verts_of_verts,
      &nedges, &verts_of_edges);
  free(verts_of_verts_offsets);
  free(verts_of_verts);
  unsigned* edges_of_verts_offsets;
  unsigned* edges_of_verts;
  up_from_down(1, 0, nedges, nverts, verts_of_edges,
      &edges_of_verts_offsets, &edges_of_verts, 0);
  unsigned* edges_of_tris = reflect_down(2, 1, ntris, nedges,
      verts_of_tris, edges_of_verts_offsets, edges_of_verts);
  free(edges_of_verts_offsets);
  free(edges_of_verts);
  unsigned* tris_of_edges_offsets;
  unsigned* tris_of_edges;
  unsigned* tris_of_edges_directions;
  up_from_down(2, 1, ntris, nedges, edges_of_tris,
      &tris_of_edges_offsets, &tris_of_edges, &tris_of_edges_directions);
  unsigned* candidates = collect_keys(2, 1, nedges,
      tris_of_edges_offsets, tris_of_edges, tris_of_edges_directions,
      bad_tris, key_of_tris);
  unsigned* elems_of_edges_offsets = tris_of_edges_offsets;
  unsigned* elems_of_edges = tris_of_edges;
  unsigned* elems_of_edges_directions = tris_of_edges_directions;
  unsigned* edges_of_elems = edges_of_tris;
  unsigned* edges_of_edges_offsets;
  unsigned* edges_of_edges;
  get_star(1, 2, nedges, elems_of_edges_offsets, elems_of_edges,
      edges_of_elems, &edges_of_edges_offsets, &edges_of_edges);
  refine_common(elem_dim, p_nelems, p_nverts, p_verts_of_elems, p_coords,
      1, nedges, verts_of_edges, candidates, edges_of_elems,
      elems_of_edges_offsets, elems_of_edges, elems_of_edges_directions,
      edges_of_edges_offsets, edges_of_edges);
  free(edges_of_elems);
  free(elems_of_edges_offsets);
  free(elems_of_edges);
  free(elems_of_edges_directions);
  free(edges_of_edges_offsets);
  free(edges_of_edges);
  free(verts_of_edges);
  free(candidates);
  return 1;
}
