#ifndef JACOBIAN_H
#define JACOBIAN_H

/* we define the element Jacobian such that
      J \xi + x_0 = x
   where J is the Jacobian, \xi is the column
   vector of parent or reference space coordinates,
   x_0 is the spatial location of the first vertex,
   and x is the resulting spacial location corresponding
   to \xi.
   This means the Jacobian ends up being:
     J = [ x_1 - x_0 | x_2 - x_0 | x_3 - x_0 ]
 */

static inline void element_jacobian(unsigned nv, double x[][3],
    double jac[3][3])
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
