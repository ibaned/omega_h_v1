#ifndef COLLAPSES_TO_VERTS_HPP
#define COLLAPSES_TO_VERTS_HPP

namespace omega_h {

struct mesh;

void valid_collapses_to_verts(struct mesh* m);
unsigned* collapsing_vertex_destinations(struct mesh* m);

}

#endif
