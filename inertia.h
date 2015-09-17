#ifndef INERTIA_H
#define INERTIA_H

void inertial_contribution(double m, double x[3], double c[3], double ic[3][3]);

unsigned eigenvals_3x3(double A[3][3], double l[]);

double inf_norm_3x3(double A[3][3]);

void eigenvector_3x3(double A[3][3], double l, double v[3]);

void least_inertial_axis(double IC[3][3], double a[3]);

#endif
