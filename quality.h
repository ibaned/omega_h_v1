#ifndef QUALITY_H
#define QUALITY_H

enum sliver_type {
  NOT_SLIVER,
  VERT_VERT_SLIVER,
  VERT_EDGE_SLIVER,
  EDGE_EDGE_SLIVER,
  VERT_FACE_SLIVER
};

typedef double (*quality_function)(double coords[][3]);

double triangle_quality(double coords[3][3]);
double tet_quality(double coords[4][3]);

extern quality_function const the_quality_functions[4];
extern quality_function const the_equal_order_quality_functions[4];

enum sliver_type triangle_sliver_type(
    double coords[3][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_edge_out);

enum sliver_type tet_sliver_type(
    double coords[4][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_ent_out);

typedef enum sliver_type (*sliver_type_function)(
    double coords[][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_out);

extern sliver_type_function const the_sliver_type_functions[4];

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

double* mesh_qualities(struct mesh* m);
double mesh_min_quality(struct mesh* m);

#endif
