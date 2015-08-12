#ifndef BAD_ELEM_KEYS
#define BAD_ELEM_KEYS

#include "quality.h"

void bad_elem_keys(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords,
    enum quality_type target,
    double qual_floor,
    double edge_ratio_floor,
    unsigned** bad_elems_out,
    unsigned** key_of_elems_out);

#endif
