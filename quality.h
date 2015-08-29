#ifndef QUALITY_H
#define QUALITY_H

enum quality_type {
  GOOD_ELEM,
  SHORT_EDGE_ELEM,
  SLIVER_ELEM,
  CAP_TET
};

typedef double (*quality_function)(double coords[][3]);

double triangle_quality(double coords[3][3]);
double tet_quality(double coords[4][3]);

extern quality_function const the_quality_functions[4];
extern quality_function const the_equal_order_quality_functions[4];

enum quality_type triangle_quality_type(
    double coords[3][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_edge_out);

enum quality_type tet_quality_type(
    double coords[4][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_ent_out);

typedef enum quality_type (*quality_type_function)(
    double coords[][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_out);

extern quality_type_function const the_quality_type_functions[4];

double* element_qualities(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords);

double min_element_quality(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords);

struct mesh;

double mesh_min_quality(struct mesh* m);

#endif
