#ifndef REFINE_H
#define REFINE_H

struct mesh;

unsigned refine_by_size(struct mesh** p_m, double qual_floor);

void uniformly_refine(struct mesh** p_m);

#endif
