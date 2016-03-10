#ifndef EVAL_FIELD_H
#define EVAL_FIELD_H

#include "include/omega_h.h"

double* eval_field(
    unsigned nents,
    double const* coords,
    unsigned ncomps,
    void (*fun)(double const* x, double* out));

struct mesh;

void mesh_eval_field(struct mesh* m, unsigned ent_dim, char const* name,
    unsigned ncomps, void (*fun)(double const* x, double* out));
void mesh_eval_field2(struct mesh* m, unsigned ent_dim, char const* name,
    unsigned ncomps, enum osh_transfer tt,
    void (*fun)(double const* x, double* out));

#endif
