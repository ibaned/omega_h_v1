#include "size.h"

#include "algebra.h"
#include "doubles.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"
#include "tag.h"

double* identity_size_field(
    unsigned nverts,
    unsigned const* vert_of_verts_offsets,
    unsigned const* vert_of_verts,
    double const* coords)
{
  double* out = LOOP_MALLOC(double, nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned first_use = vert_of_verts_offsets[i];
    unsigned end_use = vert_of_verts_offsets[i + 1];
    double edge_x[2][3];
    copy_vector(coords + i * 3, edge_x[0], 3);
    double max = 0;
    for (unsigned j = first_use; j < end_use; ++j) {
      unsigned ov = vert_of_verts[j];
      copy_vector(coords + ov * 3, edge_x[1], 3);
      double l = edge_length(edge_x);
      if (l > max)
        max = l;
    }
    out[i] = max;
  }
  return out;
}

double* element_sizes(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords)
{
  double* out = LOOP_MALLOC(double, nelems);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    double elem_coords[4][3];
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned vert = verts_of_elem[j];
      copy_vector(coords + vert * 3, elem_coords[j], 3);
    }
    out[i] = measure_entity(elem_dim, elem_coords);
  }
  return out;
}

/* TODO: create a proper cache for coordinates,
   element sizes, and element qualities, with
   invalidation when coordinates change */

double* mesh_element_sizes(struct mesh* m)
{
  return element_sizes(mesh_dim(m), mesh_count(m, mesh_dim(m)),
        mesh_ask_down(m, mesh_dim(m), 0),
        mesh_find_tag(m, 0, "coordinates")->d.f64);
}

double mesh_domain_size(struct mesh* m)
{
  double* sizes = mesh_element_sizes(m);
  double domsize = doubles_sum(sizes, mesh_count(m, mesh_dim(m)));
  loop_free(sizes);
  return domsize;
}

LOOP_KERNEL(measure_edge,
    unsigned const* verts_of_edges,
    double const* coords,
    double const* size,
    double* out)
  unsigned const* edge_vert = verts_of_edges + i * 2;
  double edge_coord[2][3];
  copy_vector(coords + edge_vert[0] * 3, edge_coord[0], 3);
  copy_vector(coords + edge_vert[1] * 3, edge_coord[1], 3);
  double length = edge_length(edge_coord);
  double desired_length = (size[edge_vert[0]] + size[edge_vert[1]]) / 2;
  out[i] = length / desired_length;
}

static double* measure_edges(
    unsigned nedges,
    unsigned const* verts_of_edges,
    double const* coords,
    double const* size)
{
  double* out = LOOP_MALLOC(double, nedges);
  LOOP_EXEC(measure_edge, nedges,
      verts_of_edges, coords, size, out);
  return out;
}

double* mesh_measure_edges_for_adapt(struct mesh* m)
{
  return measure_edges(mesh_count(m, 1), mesh_ask_down(m, 1, 0),
      mesh_find_tag(m, 0, "coordinates")->d.f64,
      mesh_find_tag(m, 0, "adapt_size")->d.f64);
}
