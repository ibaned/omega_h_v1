#ifndef ADAPT_H
#define ADAPT_H

struct mesh;

unsigned mesh_adapt(struct mesh** p_m,
    double size_ratio_floor,
    double good_qual,
    unsigned nsliver_layers,
    unsigned max_ops);

#endif
