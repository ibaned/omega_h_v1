#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <math.h>

static inline void copy_vector(
    double const a[],
    double b[],
    unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    b[i] = a[i];
}

static inline void subtract_vectors(
    double const a[],
    double const b[],
    double c[],
    unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    c[i] = a[i] - b[i];
}

static inline void add_vectors(
    double const a[],
    double const b[],
    double c[],
    unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    c[i] = a[i] + b[i];
}

static inline void swap_vectors(
    double a[],
    double b[],
    unsigned n)
{
  for (unsigned i = 0; i < n; ++i) {
    double tmp = a[i];
    a[i] = b[i];
    b[i] = tmp;
  }
}

static inline void cross_product(
    double const a[],
    double const b[],
    double c[])
{
  c[0] = a[1] * b[2]  - a[2] * b[1];
  c[1] = a[2] * b[0]  - a[0] * b[2];
  c[2] = a[0] * b[1]  - a[1] * b[0];
}

static inline double dot_product(double const a[], double const b[], unsigned n)
{
  double d = 0;
  for (unsigned i = 0; i < n; ++i)
    d += a[i] * b[i];
  return d;
}

static inline double vector_norm(double const a[], unsigned n)
{
  return sqrt(dot_product(a, a, n));
}

static inline double vector_squared_distance(
    double const a[],
    double const b[],
    unsigned n)
{
  double s = 0;
  for (unsigned i = 0; i < n; ++i) {
    double dc = a[i] - b[i];
    s += dc * dc;
  }
  return s;
}

static inline double vector_distance(
    double const a[],
    double const b[],
    unsigned n)
{
  return sqrt(vector_squared_distance(a, b, n));
}

static inline double det_3x3(double m[3][3])
{
  double tmp[3];
  cross_product(m[0], m[1], tmp);
  return dot_product(m[2], tmp, 3);
}

static inline double det_2x2(double m[2][2])
{
  return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

static inline void transp_3x3(double in[3][3], double out[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    out[i][j] = in[j][i];
}

static inline void scale_3x3(double m[3][3], double s)
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    m[i][j] *= s;
}

static inline void invert_2x2(double in[2][2], double out[2][2])
{
  double d = det_2x2(in);
  out[0][0] =  in[1][1] / d;
  out[0][1] = -in[0][1] / d;
  out[1][0] = -in[1][0] / d;
  out[1][1] =  in[0][0] / d;
}

static inline void invert_3x3(double in[3][3], double out[3][3])
{
  double tmp[3][3];
  double d = det_3x3(in);
  cross_product(in[1], in[2], tmp[0]);
  cross_product(in[2], in[0], tmp[1]);
  cross_product(in[0], in[1], tmp[2]);
  scale_3x3(tmp, 1.0 / d);
  transp_3x3(tmp, out);
}

#endif
