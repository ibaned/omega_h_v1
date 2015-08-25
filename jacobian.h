#ifndef JACOBIAN_H
#define JACOBIAN_H

static inline void element_jacobian(unsigned nv, double x[][3], double jac[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    jac[i][j] = 0;
  for (unsigned i = 0; i < nv; ++i)
  for (unsigned j = 0; j < 3; ++j)
    jac[j][i] = x[i + 1][j] - x[0][j];
}

typedef void (*jacobian_inverter)(double in[3][3], double out[3][3]);

extern jacobian_inverter const the_jacobian_inverters[4];

#endif
