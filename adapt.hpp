#ifndef ADAPT_HPP
#define ADAPT_HPP

struct mesh;

unsigned mesh_adapt(struct mesh* m,
    double size_ratio_floor,
    double good_qual,
    unsigned nsliver_layers,
    unsigned max_ops);

#endif
