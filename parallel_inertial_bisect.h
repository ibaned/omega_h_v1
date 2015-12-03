#ifndef PARALLEL_INERTIAL_BISECT_H
#define PARALLEL_INERTIAL_BISECT_H

void parallel_inertial_bisect(
    unsigned* p_n,
    double** p_coords,
    double** p_masses,
    unsigned** p_orig_ranks,
    unsigned** p_orig_ids);

void recursive_inertial_bisect(
    unsigned* p_n,
    double** p_coords,
    double** p_masses,
    unsigned** p_orig_ranks,
    unsigned** p_orig_ids);

struct mesh;

void balance_mesh_inertial(struct mesh** p_m);

#endif
