#ifndef REFINE_H
#define REFINE_H

struct mesh;

unsigned refine_by_size(struct mesh* m, double qual_floor);

void uniformly_refine(struct mesh* m);

#endif
