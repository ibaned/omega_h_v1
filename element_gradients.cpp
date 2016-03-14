#include "element_gradients.hpp"

#include <cassert>
#include <cstring>

#include "algebra.hpp"
#include "jacobian.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tables.hpp"
#include "tag.hpp"

namespace omega_h {

LOOP_KERNEL(execute,
    unsigned elem_dim,
    unsigned const* verts_of_elems,
    double const* coords_of_verts,
    unsigned ncomps,
    double const* comps_of_verts,
    unsigned ncomps_out,
    double* out,
    unsigned verts_per_elem)

  unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
  double elem_coords[4][3];
  double const* elem_comps[4];
  for (unsigned j = 0; j < verts_per_elem; ++j) {
    unsigned vert = verts_of_elem[j];
    elem_comps[j] = comps_of_verts + vert * ncomps;
    copy_vector(coords_of_verts + vert * 3, elem_coords[j], 3);
  }
  double jac[3][3];
  element_jacobian(elem_dim, elem_coords, jac);
  double jaci[3][3];
  invert_jacobian(elem_dim, jac, jaci);
  double* grad = out + i * ncomps_out;
  for (unsigned j = 0; j < ncomps_out; ++j)
    grad[j] = 0;
  for (unsigned j = 0; j < elem_dim; ++j)
  for (unsigned k = 0; k < ncomps; ++k)
  for (unsigned l = 0; l < 3; ++l)
    grad[k * 3 + l] += jaci[j][l] *
      (elem_comps[j + 1][k] - elem_comps[0][k]);
}


double* element_gradients(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords_of_verts,
    unsigned ncomps,
    double const* comps_of_verts)
{
  assert(elem_dim > 0);
  assert(ncomps > 0);
  unsigned ncomps_out = ncomps * 3;
  double* out = LOOP_MALLOC(double, ncomps_out * nelems);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  assert(verts_per_elem == elem_dim + 1);
  assert(verts_per_elem > 1);
  LOOP_EXEC(execute,nelems,
    elem_dim,
    verts_of_elems,
    coords_of_verts,
    ncomps,
    comps_of_verts,
    ncomps_out,
    out,
    verts_per_elem);
  return out;
}

struct const_tag* mesh_element_gradients(
    struct mesh* m, char const* name)
{
  struct const_tag* t = mesh_find_tag(m, 0, name);
  double* data = element_gradients(mesh_dim(m), mesh_count(m, mesh_dim(m)),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_tag(m, 0, "coordinates")->d.f64,
      t->ncomps, t->d.f64);
  static char const* prefix = "grad_";
  char* grad_name = LOOP_HOST_MALLOC(char, strlen(t->name) + strlen(prefix) + 1);
  strcpy(grad_name, prefix);
  strcat(grad_name, t->name);
  struct const_tag* out = mesh_add_tag(
      m, mesh_dim(m), TAG_F64, grad_name, t->ncomps * 3, data);
  loop_host_free(grad_name);
  return out;
}

}
