#ifndef INERTIA_H
#define INERTIA_H

void inertia_contribution(
    double m,
    double const* x,
    double const* c,
    double ic[3][3]);

void least_inertial_axis(double IC[3][3], double* a);

unsigned* local_inertia_mark(
    unsigned n,
    double const* coords,
    double const* masses);

#endif
