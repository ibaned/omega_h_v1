#ifndef REFINE_QUALITIES_H
#define REFINE_QUALITIES_H

struct mesh;

double* mesh_refine_qualities(struct mesh* m, unsigned src_dim,
    unsigned** p_candidates, double qual_floor, unsigned require_better);

#endif
