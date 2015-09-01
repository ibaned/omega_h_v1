#ifndef ADAPT_H
#define ADAPT_H

struct mesh;

void mesh_adapt(struct mesh** p_m,
    double size_ratio_floor,
    double qual_floor);

#endif
