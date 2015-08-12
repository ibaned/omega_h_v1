#ifndef QUALITY_H
#define QUALITY_H

enum tri_qual {
  GOOD_TRI,
  SHORT_EDGE_TRI,
  SLIVER_TRI
};

typedef double (*quality_function)(double coords[][3]);

double triangle_quality(double coords[3][3]);
double tet_quality(double coords[4][3]);

extern quality_function const the_quality_functions[4];

enum tri_qual triangle_quality_type(
    double coords[3][3],
    double qual_floor,
    double edge_ratio_floor,
    unsigned* key_edge_out);

#endif
