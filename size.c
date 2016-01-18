#include "size.h"

#include "algebra.h"
#include "doubles.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"
#include "tag.h"

double edge_length(double coords[2][3])
{
  return vector_distance(coords[1], coords[0], 3);
}

double triangle_area(double coords[3][3])
{
  double v[2][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return vector_norm(x, 3) / 2.0;
}

double triangle_z_area(double coords[3][3])
{
  double v[2][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return x[2] / 2.0;
}

double tet_volume(double coords[4][3])
{
  double v[3][3];
  subtract_vectors(coords[1], coords[0], v[0], 3);
  subtract_vectors(coords[2], coords[0], v[1], 3);
  subtract_vectors(coords[3], coords[0], v[2], 3);
  double x[3];
  cross_product(v[0], v[1], x);
  return dot_product(x, v[2], 3) / 6.0;
}

element_measure const the_element_measures[4] = {
  0,
  edge_length,
  triangle_area,
  tet_volume
};

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
  element_measure em = the_element_measures[elem_dim];
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    double elem_coords[4][3];
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned vert = verts_of_elem[j];
      copy_vector(coords + vert * 3, elem_coords[j], 3);
    }
    out[i] = em(elem_coords);
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
