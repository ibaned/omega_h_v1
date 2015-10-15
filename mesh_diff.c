#include "mesh_diff.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "doubles.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"
#include "tag.h"

static unsigned doubles_diff(double const* a, double const* b, unsigned n,
    double tol, double floor)
{
  assert(0 <= tol);
  assert(0 <= floor);
  double* diffs = LOOP_MALLOC(double, n);
  for (unsigned i = 0; i < n; ++i) {
    double fa = fabs(a[i]);
    double fb = fabs(b[i]);
    if (fa < floor && fb < floor) {
      diffs[i] = 0;
      continue;
    }
    double fm = fb > fa ? fb : fa;
    diffs[i] = fabs(a[i] - b[i]) / fm;
  }
  double maxdiff = doubles_max(diffs, n);
  loop_free(diffs);
  if (maxdiff > tol) {
    printf("max relative difference %e\n", maxdiff);
    return 1;
  }
  return 0;
}

static unsigned uchars_diff(unsigned char const* a, unsigned char const* b, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] != b[i])
      return 1;
  return 0;
}

static unsigned uints_diff(unsigned const* a, unsigned const* b, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] != b[i])
      return 1;
  return 0;
}

static unsigned ulongs_diff(unsigned long const* a, unsigned long const* b, unsigned n)
{
  for (unsigned i = 0; i < n; ++i)
    if (a[i] != b[i])
      return 1;
  return 0;
}

static unsigned tag_diff(struct const_tag* a, struct const_tag* b, unsigned n,
    double tol, double floor)
{
  if (strcmp(a->name, b->name)) {
    printf("tag names \"%s\" and \"%s\" differ\n",
        a->name, b->name);
    return 1;
  }
  if (a->ncomps != b->ncomps) {
    printf("tag \"%s\" sizes %u and %u differ\n",
        a->name, a->ncomps, b->ncomps);
    return 1;
  }
  if (a->type != b->type) {
    printf("tag \"%s\" types differ\n", a->name);
    return 1;
  }
  switch (a->type) {
    case TAG_U8:
      if (uchars_diff(a->d.u8, b->d.u8, n * a->ncomps)) {
        printf("tag \"%s\" contents differ\n", a->name);
        return 1;
      }
      break;
    case TAG_U32:
      if (uints_diff(a->d.u32, b->d.u32, n * a->ncomps)) {
        printf("tag \"%s\" contents differ\n", a->name);
        return 1;
      }
      break;
    case TAG_U64:
      if (ulongs_diff(a->d.u64, b->d.u64, n * a->ncomps)) {
        printf("tag \"%s\" contents differ\n", a->name);
        return 1;
      }
      break;
    case TAG_F64:
      if (doubles_diff(a->d.f64, b->d.f64, n * a->ncomps, tol, floor)) {
        printf("tag \"%s\" contents differ\n", a->name);
        return 1;
      }
      break;
  }
  return 0;
}

static unsigned tags_diff(struct tags* a, struct tags* b, unsigned n,
    double tol, double floor)
{
  unsigned ntags = b->n > a->n ? b->n : a->n;
  for (unsigned i = 0; i < ntags; ++i)
    if (tag_diff(get_tag(a, i), get_tag(b, i), n, tol, floor))
      return 1;
  if (a->n == b->n)
    return 0;
  for (unsigned i = ntags; i < a->n; ++i)
    printf("tag \"%s\" only on first mesh\n",
        get_tag(a, i)->name);
  for (unsigned i = ntags; i < b->n; ++i)
    printf("tag \"%s\" only on second mesh\n",
        get_tag(b, i)->name);
  return 1;
}

unsigned mesh_diff(struct mesh* a, struct mesh* b,
    double tol, double floor)
{
  if (mesh_dim(a) != mesh_dim(b)) {
    printf("mesh dimensionalities %u and %u differ\n",
        mesh_dim(a), mesh_dim(b));
    return 1;
  }
  unsigned dim = mesh_dim(a);
  for (unsigned d = 0; d <= dim; ++d) {
    if (!mesh_has_dim(a, d) && !mesh_has_dim(b, d))
      continue;
    if (mesh_has_dim(a, d) && !mesh_has_dim(b, d)) {
      printf("only first mesh has dimension %u\n", d);
      return 1;
    }
    if (!mesh_has_dim(a, d) && mesh_has_dim(b, d)) {
      printf("only second mesh has dimension %u\n", d);
      return 1;
    }
    if (mesh_count(a, d) != mesh_count(b, d)) {
      printf("counts %u and %u differ for dimension %u\n",
          mesh_count(a, d), mesh_count(b, d), d);
      return 1;
    }
    unsigned nents = mesh_count(a, d);
    if (d) {
      unsigned const* eva = mesh_ask_down(a, d, 0);
      unsigned const* evb = mesh_ask_down(b, d, 0);
      unsigned vpe = the_down_degrees[d][0];
      if (uints_diff(eva, evb, vpe * nents)) {
        printf("connectivity differs for dimension %u\n", d);
        return 1;
      }
    }
    struct tags* tsa = mesh_tags(a, d);
    struct tags* tsb = mesh_tags(b, d);
    if (tags_diff(tsa, tsb, nents, tol, floor)) {
      printf("tags for dimension %u differ\n", d);
      return 1;
    }
  }
  return 0;
}
