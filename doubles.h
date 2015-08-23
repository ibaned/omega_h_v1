#ifndef DOUBLES_H
#define DOUBLES_H

double doubles_max(double const a[], unsigned n);
double doubles_min(double const a[], unsigned n);
void doubles_axpy(double a, double const x[], double const y[],
    double out[], unsigned n);

#endif
