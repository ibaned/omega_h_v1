#ifndef FIND_ROOTS_H
#define FIND_ROOTS_H

/* finds the real roots of the linear polynomial
   a * x + b = 0
   the return value is the number of roots */
unsigned find_linear_root(double a, double b, double* root);

/* finds the real roots of the quadratic polynomial
   a * x^2 + b * x + c = 0
   the return value is the number of roots */
unsigned find_quadratic_roots(
    double a, double b, double c,
    double* roots);

/* finds the real roots of the cubic polynomial
   a * x^3 + b * x^2 + c * x + d = 0
   the return value is the number of roots */
unsigned find_cubic_roots(
    double a, double b, double c, double d,
    double* roots);

#endif
