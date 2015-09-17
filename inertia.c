#include "inertia.h"
#include "algebra.h"
#include "find_roots.h"
#include <assert.h>

#include <stdio.h>

static void cross_matrix(double b[3], double B[3][3])
{
  B[0][0] =     0; B[0][1] = -b[2]; B[0][2] =  b[1];
  B[1][0] =  b[2]; B[1][1] =     0; B[1][2] = -b[0];
  B[2][0] = -b[1]; B[2][1] =  b[0]; B[2][2] =    0;
}

static void mmm_3x3(double a[3][3], double b[3][3], double c[3][3])
{
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j) {
    c[i][j] = 0;
    for (unsigned k = 0; k < 3; ++k)
      c[i][j] += a[i][k] * b[k][j];
  }
}

void inertial_contribution(double m, double x[3], double c[3], double ic[3][3])
{
  double dx[3];
  subtract_vectors(x, c, dx, 3);
  double B[3][3];
  cross_matrix(dx, B);
  mmm_3x3(B, B, ic);
  scale_3x3(ic, -m, ic);
}

static double trace_3x3(double A[3][3])
{
  return A[0][0] + A[1][1] + A[2][2];
}

double inf_norm_3x3(double A[3][3])
{
  double maxv = 0;
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    if (fabs(A[i][j]) > maxv)
      maxv = fabs(A[i][j]);
  return maxv;
}

static void char_poly_3x3(double A[3][3],
    double* a, double* b, double* c, double* d)
{
  *a = -1;
  double tA = trace_3x3(A);
  *b = tA;
  double Asq[3][3];
  mmm_3x3(A, A, Asq);
  *c = -(tA * tA - trace_3x3(Asq)) / 2.;
  *d = det_3x3(A);
}

unsigned eigenvals_3x3(double A[3][3], double l[])
{
  double a,b,c,d;
  char_poly_3x3(A, &a, &b, &c, &d);
  return find_cubic_roots(a, b, c, d, l);
}

/* returns a vector (o) orthogonal to vector (v).
   this function should be robust to any value of (v) */
static void ortho_vector(double v[3], double o[3])
{
  unsigned a = 0, b = 1, c = 2;
  unsigned tmp;
#define S(i,j) if (fabs(v[i]) < fabs(v[j])) {tmp = i;i = j;j = tmp;}
  S(a,b); /* "handcoded" sort for three numbers... */
  S(b,c);
  S(a,b);
#undef S
  /* now v[a], v[b], v[c] should decrease in magnitude */
  assert(fabs(v[a]) >= fabs(v[b]));
  assert(fabs(v[b]) >= fabs(v[c]));
  o[a] = - v[b]; /* 90 degree rotation in a-b plane */
  o[b] = v[a];
  o[c] = 0; /* zero component along c */
  /* note: we may want to scale this such that
     the norm of (o) is the same as the norm of (v),
     right now the loss of the v[c] term breaks that */
}

/* finds an eigenvector of matrix (A) associated with eigenvalue (l) */
void eigenvector_3x3(double A[3][3], double l, double v[3])
{
  double s[3][3];
  transp_3x3(A, s);
  for (unsigned i = 0; i < 3; ++i)
    s[i][i] -= l;
  /* s = transpose(A) - l * I */
  /* s is assumed to be rank-deficient.
     it may have rank 2, 1, or zero */
  assert(det_3x3(s) < 1e-10);
  /* we are looking for a vector v such that
     s * v = 0
     which means the v must be orthogonal to
     every row of s, hence the transpose */
  /* we first assume it has rank 2.
     we may be able to capture the null space
     vector by taking the largest cross product
     of the rows */
  unsigned cpi[3][2] = {{0,1},{1,2},{0,2}};
  v[0]=v[1]=v[2]=0;
  for (unsigned i = 0; i < 3; ++i) {
    double cp[3];
    cross_product(s[cpi[i][0]],s[cpi[i][1]],cp);
    if (vector_norm(cp,3) > vector_norm(v,3))
      copy_vector(cp, v, 3);
  }
  if (vector_norm(v,3) > 1e-10) { /*these checks can be argued about*/
    printf("rank 2\n");
    return;
  }
  /* okay, maybe not rank 2. maybe rank 1.
     the null space is a plane orthogonal
     to all rows of s, so lets take the largest
     row of s and try a vector orthogonal to that */
  for (unsigned i = 0; i < 3; ++i)
    if (vector_norm(s[i],3) > vector_norm(v,3))
      ortho_vector(s[i], v);
  if (vector_norm(v,3) > 1e-10) { /*these checks can be argued about*/
    printf("rank 1\n");
    return;
  }
  /* okay, this is a zero matrix (rank 0).
     the null space is everywhere.
     pick any vector, in particular the most obvious one */
  v[0] = 1;
  v[1] = 0;
  v[2] = 0;
  printf("rank 0\n");
}
