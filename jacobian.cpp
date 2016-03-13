#include "jacobian.hpp"

#include "algebra.hpp"

namespace omega_h {

LOOP_IN void
element_jacobian(unsigned nv, double (*x)[3],
    double jac[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    jac[i][j] = 0;
  for (unsigned i = 0; i < nv; ++i)
  for (unsigned j = 0; j < 3; ++j)
    jac[j][i] = x[i + 1][j] - x[0][j];
}

LOOP_IN void
invert_jacobian_1(double in[3][3], double out[3][3])
{
  double msq = dot_product(in[0], in[0], 3);
  for (unsigned i = 0; i < 3; ++i)
    out[i][0] = in[0][i] / msq;
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 1; j < 3; ++j)
    out[i][j] = 0;
}

LOOP_IN void
invert_jacobian_2(double in[3][3], double out[3][3])
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

LOOP_IN void
invert_jacobian_3(double in[3][3], double out[3][3])
{
  invert_3x3(in, out);
}

LOOP_IN void
invert_jacobian(unsigned dim, double in[3][3], double out[3][3])
{
  switch (dim) {
    case 3: invert_jacobian_3(in, out); return;
    case 2: invert_jacobian_2(in, out); return;
    case 1: invert_jacobian_1(in, out); return;
  }
  LOOP_NORETURN();
}

}
