#ifndef SLIVER_KEYS_H
#define SLIVER_KEYS_H

#include "quality.h"

void sliver_keys(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords,
    enum sliver_type target,
    double qual_floor,
    double edge_ratio_floor,
    unsigned** bad_elems_out,
    unsigned** key_of_elems_out);

#endif
