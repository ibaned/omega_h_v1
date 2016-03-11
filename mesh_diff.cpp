#include "mesh_diff.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "doubles.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tables.hpp"
#include "tag.hpp"

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
    for (unsigned i = 0; i < n; ++i) {
      if (diffs[i] == maxdiff) {
        printf("max difference pair is (%e, %e), #%u\n",
            a[i], b[i], i);
      }
    }
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
    double tol, double floor, unsigned allow_superset)
{
  for (unsigned i = 0; i < count_tags(a); ++i) {
    struct const_tag* ta = get_tag(a, i);
    struct const_tag* tb = find_tag(b, ta->name);
    if (!tb) {
      if (allow_superset)
        continue;
      printf("tag \"%s\" only on first mesh\n", ta->name);
      return 1;
    }
    if (tag_diff(ta, tb, n, tol, floor))
      return 1;
  }
  if (!allow_superset)
    for (unsigned i = 0; i < count_tags(b); ++i)
      if (!find_tag(a, get_tag(b, i)->name)) {
        printf("tag \"%s\" only on second mesh\n", get_tag(b, i)->name);
        return 1;
      }
  return 0;
}

unsigned mesh_diff(struct mesh* a, struct mesh* b,
    double tol, double floor, unsigned allow_superset)
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
      if (allow_superset)
        continue;
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
    if (tags_diff(tsa, tsb, nents, tol, floor, allow_superset)) {
      printf("tags for dimension %u differ\n", d);
      return 1;
    }
  }
  return 0;
}
