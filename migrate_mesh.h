#ifndef MIGRATE_MESH_H
#define MIGRATE_MESH_H

struct mesh;

void migrate_mesh(struct mesh** p_m,
    unsigned nelems_recvd,
    unsigned const* recvd_elem_ranks,
    unsigned const* recvd_elem_ids);

#endif
