#ifndef SIZE_H
#define SIZE_H

typedef double (*size_function)(double const coords[][3]);

double edge_length(double const coords[2][3]);
double triangle_area(double const coords[3][3]);
double tet_volume(double const coords[4][3]);

extern size_function const the_size_functions[4];

#endif
