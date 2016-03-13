#ifndef REFINE_NODAL_HPP
#define REFINE_NODAL_HPP

namespace omega_h {

double* refine_nodal(
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* gen_offset_of_srcs,
    unsigned comps_per_vert,
    double const* field);

}

#endif
