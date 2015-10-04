#ifndef INERTIA_H
#define INERTIA_H

void inertial_contribution(double m, double* x, double* c, double ic[3][3]);

void least_inertial_axis(double IC[3][3], double* a);

#endif
