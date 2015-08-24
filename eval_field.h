#ifndef EVAL_FIELD_H
#define EVAL_FIELD_H

double* eval_field(
    unsigned nverts,
    double const* coords,
    unsigned ncomps,
    void (*fun)(double const x[3], double out[]));

struct mesh;

void mesh_eval_field(struct mesh* m, char const* name, unsigned ncomps,
    void (*fun)(double const x[3], double out[]));

#endif
