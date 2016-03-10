#ifndef COLLAPSES_TO_VERTS_H
#define COLLAPSES_TO_VERTS_H

struct mesh;

void valid_collapses_to_verts(struct mesh* m);
unsigned* collapsing_vertex_destinations(struct mesh* m);

#endif
