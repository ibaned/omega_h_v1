#ifndef INERTIA_H
#define INERTIA_H

void inertial_contribution(double m, double x[3], double c[3], double ic[3][3]);

void least_inertial_axis(double ic[3][3], double a[]);

#endif
