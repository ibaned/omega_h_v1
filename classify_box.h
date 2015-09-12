#ifndef CLASSIFY_BOX_H
#define CLASSIFY_BOX_H

struct mesh;

unsigned* classify_box(
    unsigned nverts,
    double const* coords);

void mesh_classify_box(struct mesh* m);

#endif
