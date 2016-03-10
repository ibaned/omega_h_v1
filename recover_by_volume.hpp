#ifndef RECOVER_BY_VOLUME_H
#define RECOVER_BY_VOLUME_H

double* recover_by_volume(
    unsigned nverts,
    unsigned const* elems_of_verts_offsets,
    unsigned const* elems_of_verts,
    double const* size_of_elems,
    unsigned ncomps,
    double const* comps_of_elems);

struct mesh;

struct const_tag* mesh_recover_by_volume(
    struct mesh* m, char const* name);

#endif
