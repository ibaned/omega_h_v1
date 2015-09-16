#include "find_roots.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* kicks off the insert_root process with
   a[0] * x + a[1] = (x - r)        */
static void init_root(double a[], double r)
{
  a[0] = 1;
  a[1] = -r;
  printf("init root %f %f\n", a[0], a[1]);
}

/* this function multiplies the polynomial
  a[0] * x ^ (n - 1) + ... + a[n - 1] * x + a[n]
times the binomial (x - r), resulting in the polynomial
  a[0] * x ^ n + ... + a[n] * x + a[n + 1]        */
static void insert_root(double a[], unsigned n, double r)
{
  a[n] = 0;
  for (unsigned i = 0; i <= n; ++i)
    a[i + 1] -= a[i] * r;
  printf("insert root");
  for (unsigned i = 0; i <= n + 1; ++i)
    printf(" %f", a[i]);
  printf("\n");
}

/* returns a uniform random real number from 0 to 1 */
static double real_rand(void)
{
  return ((double)(rand())) / ((double)RAND_MAX);
}

/* chooses one of a set of cases based on their
   probabilities. probs[i] has the probability
   of case (i) being selected, and the sum
   of probabilities should be 1 */
static unsigned case_rand(double probs[], unsigned ncases)
{
  double v = real_rand();
  for (unsigned i = 0; i < ncases; ++i) {
    if (v < probs[i])
      return i;
    v -= probs[i];
  }
  return ncases - 1;
}

static unsigned gen_test(double a[], double r[])
{
  /* probabilities of testing a given
     polynomial order */
  double dim_probs[4] = {0.0,0.1,0.2,0.7};
  unsigned dim = case_rand(dim_probs, 4);
  printf("picking dimension %u\n", dim);
  printf("known roots");
  for (unsigned i = 0; i < dim; ++i) {
    /* probability of a zero root */
    if (real_rand() < 0.1)
      r[i] = 0;
    else
      r[i] = real_rand();
    printf(" %f", r[i]);
  }
  printf("\n");
  init_root(a, r[0]);
  for (unsigned i = 1; i < dim; ++i)
    insert_root(a, i, r[i]);
  return dim;
}

static void run_test(void)
{
  double a[4];
  double r[3];
  unsigned dim = gen_test(a, r);
  double fr[3];
  for (unsigned i = 0; i <= dim; ++i)
    a[i + (3 - dim)] = a[i];
  for (unsigned i = 0; i < (3 - dim); ++i)
    a[i] = 0;
  printf("coefficients:");
  for (unsigned i = 0; i < 4; ++i)
    printf(" %f", a[i]);
  printf("\n");
  unsigned nroots = find_cubic_roots(
      a[0], a[1], a[2], a[3], fr);
  printf("found roots");
  for (unsigned i = 0; i < nroots; ++i)
    printf(" %f", fr[i]);
  printf("\n");
  for (unsigned i = 0; i < dim; ++i) {
    double known_root = r[i];
    unsigned j;
    for (j = 0; j < nroots; ++j) {
      double found_root = fr[j];
      if (fabs(found_root - known_root) < 1e-10)
        break;
    }
    if (j == nroots)
      abort();
  }
}

int main()
{
  run_test();
}
