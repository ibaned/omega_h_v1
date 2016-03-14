#ifndef COARSEN_TOPOLOGY_HPP
#define COARSEN_TOPOLOGY_HPP

namespace omega_h {

void coarsen_topology(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* gen_offset_of_elems,
    unsigned const* gen_vert_of_elems,
    unsigned const* gen_direction_of_elems,
    unsigned* ngen_elems_out,
    unsigned** verts_of_gen_elems_out);

}

#endif
