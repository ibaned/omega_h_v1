#ifndef ELEMENT_GRADIENTS_HPP
#define ELEMENT_GRADIENTS_HPP

namespace omega_h {

double* element_gradients(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords_of_verts,
    unsigned ncomps,
    double const* comps_of_verts);

struct mesh;

struct const_tag* mesh_element_gradients(
    struct mesh* m, char const* name);

}

#endif
