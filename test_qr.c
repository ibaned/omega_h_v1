#include "qr.h"
#include "algebra.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

#define TOLERANCE 1e-10

static double ident[3][3] = {
  {1,0,0},
  {0,1,0},
  {0,0,1}
};

static void print(double a[3][3], char const* name)
{
  printf("%s:\n", name);
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j)
      printf(" %f", a[i][j]);
    printf("\n");
  }
}

static void assert_close(double a[3][3], double b[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    assert(fabs(a[i][j] - b[i][j]) < TOLERANCE);
}

static void verify_q(double q[3][3])
{
  double qt[3][3];
  transp_3x3(q,qt);
  print(qt, "Q^T");
  double qqt[3][3];
  mul_3x3(q,qt,qqt);
  print(qqt, "Q*Q^T");
  assert_close(qqt, ident);
}

static void verify_r(double r[3][3])
{
  for (unsigned j = 0; j < 3; ++j)
  for (unsigned i = j + 1; i < 3; ++i)
    assert(fabs(r[i][j]) < TOLERANCE);
}

static void test_one_qr(double a[3][3])
{
  print(a, "A"); 
  double q[3][3];
  double r[3][3];
  qr_decomp(a,q,r);
  print(q, "Q");
  print(r, "R");
  double qr[3][3];
  mul_3x3(q,r,qr);
  print(qr, "QR");
  assert_close(qr, a);
  verify_q(q);
  verify_r(r);
}

int main()
{
  double a[3][3] = {
    {0,0,0},
    {0,1,0},
    {0,0,0}
  };
  test_one_qr(a);
}
