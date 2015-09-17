#ifndef INERTIA_H
#define INERTIA_H

void inertial_contribution(double m, double x[3], double c[3], double ic[3][3]);

unsigned eigenvals_3x3(double A[3][3], double l[]);

#endif
