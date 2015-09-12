#include "element_gradients.h"
#include <string.h>    // for strcat, strcpy
#include <stdlib.h>    // for malloc, free
#include <string.h>    // for strlen
#include <assert.h>    // for assert
#include "algebra.h"   // for copy_vector
#include "field.h"     // for const_field
#include "jacobian.h"  // for element_jacobian, jacobian_inverter, the...
#include "mesh.h"      // for mesh_dim, mesh_find_nodal_field, mesh_ad...
#include "tables.h"    // for the_down_degrees

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
  double* out = malloc(sizeof(double) * ncomps_out * nelems);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  assert(verts_per_elem == elem_dim + 1);
  jacobian_inverter jif = the_jacobian_inverters[elem_dim];
  for (unsigned i = 0; i < nelems; ++i) {
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
    jif(jac, jaci);
    double* grad = out + i * ncomps_out;
    for (unsigned j = 0; j < ncomps_out; ++j)
      grad[j] = 0;
    for (unsigned j = 0; j < elem_dim; ++j)
    for (unsigned k = 0; k < ncomps; ++k)
    for (unsigned l = 0; l < 3; ++l) {
      assert(elem_comps[j + 1]);
      assert(elem_comps[0]);
      unsigned kl = k * 3 + l;
      assert(l < ncomps_out);
      grad[kl] += (elem_comps[j + 1][k] - elem_comps[0][k]) * jaci[j][l];
    }
  }
  return out;
}

struct const_field* mesh_element_gradients(
    struct mesh* m, char const* name)
{
  struct const_field* f = mesh_find_nodal_field(m, name);
  double* data = element_gradients(mesh_dim(m), mesh_count(m, mesh_dim(m)),
      mesh_ask_down(m, mesh_dim(m), 0),
      mesh_find_nodal_field(m, "coordinates")->data,
      f->ncomps, f->data);
  static char const* prefix = "grad_";
  char* grad_name = malloc(strlen(f->name) + strlen(prefix) + 1);
  strcpy(grad_name, prefix);
  strcat(grad_name, f->name);
  struct const_field* out = mesh_add_elem_field(
      m, grad_name, f->ncomps * 3, data);
  free(grad_name);
  return out;
}
