#ifndef QR_HPP
#define QR_HPP

#include "loop.hpp"

#define MAX_PTS 16

namespace omega_h {

void qr_decomp(double a[3][3], double q[3][3], double r[3][3]);
LOOP_INOUT unsigned qr_decomp2(
    double a[MAX_PTS][4],
    double q[MAX_PTS][MAX_PTS],
    double r[MAX_PTS][4],
    unsigned npts);
LOOP_INOUT void qr_solve2(
    double q[MAX_PTS][MAX_PTS],
    double r[MAX_PTS][4],
    double b[MAX_PTS],
    unsigned npts, double x[4]);
void least_squares_fit(double a[MAX_PTS][4], double b[MAX_PTS],
    unsigned npts, double x[4]);
void qr_eigen(double a[3][3], double q[3][3], double l[3][3]);

}

#endif
