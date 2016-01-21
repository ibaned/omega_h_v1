#ifndef REFINE_QUALITIES_H
#define REFINE_QUALITIES_H

double* refine_qualities(
    unsigned elem_dim,
    unsigned src_dim,
    unsigned nsrcs,
    unsigned const* verts_of_srcs,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_srcs_offsets,
    unsigned const* elems_of_srcs,
    unsigned const* elems_of_srcs_directions,
    /* candidates are rejected whose output qualities
       would be worse than (qual_floor),
       and in the case of (require_better), output qualities
       worse than input qualities */
    unsigned* candidate_srcs,
    double const* coords,
    double qual_floor,
    double const* elem_quals,
    unsigned require_better);

struct mesh;

double* mesh_refine_qualities(struct mesh* m, unsigned src_dim,
    unsigned** p_candidates, double qual_floor, unsigned require_better);

#endif
