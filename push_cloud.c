#include "push_cloud.h"
#include "cloud.h"
#include "mesh.h"
#include "tables.h"
#include "algebra.h"
#include "ints.h"
#include "loop.h"
#include <assert.h>

static void push_particle_2d(
    unsigned pt,
    unsigned* elem_of_pts,
    unsigned* last_dir_of_pts,
    double const* coords_of_pts,
    unsigned const* verts_of_elems,
    double const* coords_of_verts,
    unsigned const* dual)
{
  double const* coords_of_pt = coords_of_pts + pt * 3;
  unsigned const* const* elem_verts_of_edges =
    the_canonical_orders[2][1][0];
  unsigned elem = elem_of_pts[pt];
  unsigned dir;
  while (1) {
    unsigned const* verts_of_elem = verts_of_elems + elem * 3;
    dir = INVALID;
    for (unsigned i = 0; i < 3; ++i) {
      unsigned const* elem_verts_of_edge =
        elem_verts_of_edges[i];
      double ex[2][3];
      for (unsigned j = 0; j < 2; ++j) {
        unsigned elem_vert = elem_verts_of_edge[j];
        unsigned vert = verts_of_elem[elem_vert];
        copy_vector(coords_of_verts + vert * 3, ex[j], 3);
      }
      double ev[3];
      subtract_vectors(ex[1], ex[0], ev, 3);
      double nl[3];
      nl[0] = -ev[1];
      nl[1] = ev[0];
      nl[2] = 0;
      double rel[3];
      subtract_vectors(coords_of_pt, ex[0], rel, 3);
      double d = dot_product(nl, rel, 3);
      if (d < 0) {
        dir = i;
        break;
      }
    }
    if (dir == INVALID) {
      break;
    } else {
      unsigned next_elem = dual[elem * 3 + dir];
      if (next_elem == INVALID)
        break;
      elem = next_elem;
    }
  }
  elem_of_pts[pt] = elem;
  last_dir_of_pts[pt] = dir;
}

void push_cloud(struct cloud* c, struct mesh* m)
{
  assert(mesh_dim(m) == 2);
  unsigned npts = cloud_count(c);
  unsigned* elem_of_pts = ints_copy(cloud_find_tag(c, "mesh_elem")->data, npts);
  unsigned* last_dir_of_pts = loop_malloc(sizeof(unsigned) * npts);
  double const* coords_of_pts = cloud_find_tag(c, "coordinates")->data;
  unsigned const* verts_of_elems = mesh_ask_down(m, mesh_dim(m), 0);
  double const* coords_of_verts = mesh_find_tag(m, 0, "coordinates")->data;
  unsigned const* dual = mesh_ask_dual(m);
  for (unsigned i = 0; i < npts; ++i)
    push_particle_2d(i, elem_of_pts, last_dir_of_pts, coords_of_pts,
        verts_of_elems, coords_of_verts, dual);
  cloud_free_tag(c, "mesh_elem");
  cloud_add_tag(c, TAG_U32, "mesh_elem", 1, elem_of_pts);
  cloud_add_tag(c, TAG_U32, "last_dir", 1, last_dir_of_pts);
}
