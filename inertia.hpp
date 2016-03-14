#ifndef INERTIA_HPP
#define INERTIA_HPP

namespace omega_h {

void inertia_contribution(
    double m,
    double const* x,
    double const* c,
    double ic[3][3]);

void least_inertial_axis(double IC[3][3], double* a);

unsigned* mark_inertial_bisection(
    unsigned n,
    double const* coords,
    double const* masses,
    unsigned is_global);

}

#endif
