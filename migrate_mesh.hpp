#ifndef MIGRATE_MESH_HPP
#define MIGRATE_MESH_HPP

namespace omega_h {

struct mesh;
struct exchanger;

void migrate_mesh(struct mesh* m,
    unsigned nelems_recvd,
    unsigned const* recvd_elem_ranks,
    unsigned const* recvd_elem_ids);

void migrate_mesh(struct mesh* m, struct exchanger* elem_push);

}

#endif
