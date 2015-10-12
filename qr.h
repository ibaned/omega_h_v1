#ifndef QR_H
#define QR_H

void qr_decomp(double a[3][3], double q[3][3], double r[3][3]);
void qr_eigen(double a[3][3], double q[3][3], double l[3][3]);

#endif