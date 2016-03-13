#ifndef JACOBIAN_HPP
#define JACOBIAN_HPP

#include "algebra.hpp"

namespace omega_h {

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

LOOP_IN void
element_jacobian(unsigned nv, double (*x)[3],
    double jac[3][3]);

LOOP_IN void
invert_jacobian(unsigned dim, double in[3][3], double out[3][3]);

}

#endif
