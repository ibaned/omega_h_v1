#include "quality.hpp"

#include "arrays.hpp"
#include "doubles.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "size.hpp"
#include "tables.hpp"
#include "tag.hpp"

LOOP_KERNEL(elem_quality_kern,
    unsigned const* verts_of_elems,
    unsigned elem_dim,
    unsigned verts_per_elem,
    double const* coords,
    double* out)
  unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
  double elem_x[MAX_DOWN][3];
  for (unsigned j = 0; j < verts_per_elem; ++j) {
    unsigned vert = verts_of_elem[j];
    copy_vector(coords + vert * 3, elem_x[j], 3);
  }
  out[i] = element_quality(elem_dim, elem_x);
}

double* element_qualities(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords)
{
  if (elem_dim < 2)
    return filled_array(nelems, 1.0);
  double* out = LOOP_MALLOC(double, nelems);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  LOOP_EXEC(elem_quality_kern, nelems,
      verts_of_elems,
      elem_dim,
      verts_per_elem,
      coords,
      out);
  return out;
}

double min_element_quality(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords)
{
  double* quals = element_qualities(elem_dim, nelems, verts_of_elems, coords);
  double mq = doubles_min(quals, nelems);
  loop_free(quals);
  return mq;
}

double* mesh_qualities(struct mesh* m)
{
  return element_qualities(mesh_dim(m),
      mesh_count(m, mesh_dim(m)),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_tag(m, 0, "coordinates")->d.f64);
}

double mesh_min_quality(struct mesh* m)
{
  return min_element_quality(mesh_dim(m),
      mesh_count(m, mesh_dim(m)),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_tag(m, 0, "coordinates")->d.f64);
}
