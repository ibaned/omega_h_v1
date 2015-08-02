#ifndef QUALITY_H
#define QUALITY_H

typedef double (*quality_function)(double const coords[][3]);

double triangle_quality(double const coords[3][3]);
double tet_quality(double const coords[4][3]);

extern quality_function const the_quality_functions[4];

#endif
