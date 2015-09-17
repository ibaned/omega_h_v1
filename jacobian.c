#include "jacobian.h"
#include "algebra.h"

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


static void invert_jacobian_1(double in[3][3], double out[3][3])
{
  double msq = dot_product(in[0], in[0], 3);
  for (unsigned i = 0; i < 3; ++i)
    out[i][0] = in[0][i] / msq;
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 1; j < 3; ++j)
    out[i][j] = 0;
}

static void invert_jacobian_2(double in[3][3], double out[3][3])
{
  double normal[2][2];
  for (unsigned i = 0; i < 2; ++i)
  for (unsigned j = 0; j < 2; ++j) {
    normal[i][j] = 0;
    for (unsigned k = 0; k < 3; ++k)
      normal[i][j] += in[i][k] * in[j][k];
  }
  double inv_normal[2][2];
  invert_2x2(normal, inv_normal);
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    out[i][j] = 0;
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 2; ++j)
  for (unsigned k = 0; k < 2; ++k)
    out[i][j] += in[k][i] * inv_normal[k][j];
}

static void invert_jacobian_3(double in[3][3], double out[3][3])
{
  invert_3x3(in, out);
}

jacobian_inverter const the_jacobian_inverters[4] = {
  0,
  invert_jacobian_1,
  invert_jacobian_2,
  invert_jacobian_3
};
