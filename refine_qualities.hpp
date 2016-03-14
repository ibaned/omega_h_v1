#ifndef REFINE_QUALITIES_HPP
#define REFINE_QUALITIES_HPP

namespace omega_h {

struct mesh;

double* mesh_refine_qualities(struct mesh* m, unsigned src_dim,
    unsigned** p_candidates, double qual_floor, unsigned require_better);

}

#endif
