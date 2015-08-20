#ifndef INSIDE_H
#define INSIDE_H

void point_to_edge(double edge[2][3], double const pt[3], double b[2]);
void point_to_triangle(double tri[3][3], double const pt[3], double b[3]);
void point_to_tet(double tet[4][3], double const pt[3], double b[4]);

typedef void (*inside_function)(double elem[][3], double const pt[3], double b[]);

extern inside_function const the_inside_functions[4];

#endif
