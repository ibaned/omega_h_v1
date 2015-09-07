#ifndef REFINE_BY_SIZE_H
#define REFINE_BY_SIZE_H

struct mesh;

unsigned refine_by_size(struct mesh** p_m, double qual_floor);

#endif
