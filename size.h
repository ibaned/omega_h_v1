#ifndef SIZE_H
#define SIZE_H

typedef double (*element_measure)(double coords[][3]);

double edge_length(double coords[2][3]);
double triangle_area(double coords[3][3]);
double triangle_xy_area(double coords[3][3]);
double tet_volume(double coords[4][3]);

extern element_measure const the_element_measures[4];

double* size_function_to_field(
    unsigned nverts,
    double const* coords,
    double (*f)(double const x[]));

double* identity_size_field(
    unsigned nverts,
    unsigned const* vert_of_verts_offsets,
    unsigned const* vert_of_verts,
    double const* coords);

double* clamp_size_field(
    unsigned nverts,
    double const* id_size,
    double const* user_size,
    double ceil,
    double floor,
    unsigned* did_clamp_out);

#endif
