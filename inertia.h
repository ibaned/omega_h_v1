#ifndef INERTIA_H
#define INERTIA_H

void inertial_contribution(
    double m,
    double const* x,
    double const* c,
    double ic[3][3]);

void least_inertial_axis(double IC[3][3], double* a);

void local_inertial_mark(
    unsigned n,
    double const* coords,
    double const* masses,
    double tol,
    unsigned** in);

#endif
