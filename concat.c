#include "concat.h"
#include "tables.h"
#include <stdlib.h>
#include <string.h>

static void* concat_general(
    unsigned typesize,
    unsigned width,
    void const* a,
    unsigned na,
    void const* b,
    unsigned nb)
{
  unsigned a_bytes = na * width * typesize;
  unsigned b_bytes = na * width * typesize;
  unsigned total_bytes = a_bytes + b_bytes;
  char* out = malloc(total_bytes);
  memcpy(out, a, a_bytes);
  memcpy(out + a_bytes, b, b_bytes);
  return out;
}

unsigned* concat_ints(
    unsigned width,
    unsigned const a[],
    unsigned na,
    unsigned const b[],
    unsigned nb)
{
  return concat_general(sizeof(unsigned), width, a, na, b, nb);
}

double* concat_doubles(
    unsigned width,
    double const a[],
    unsigned na,
    double const b[],
    unsigned nb)
{
  return concat_general(sizeof(double), width, a, na, b, nb);
}

void* general_subset(
    unsigned typesize,
    unsigned n,
    unsigned width,
    void const* a,
    unsigned const* offsets)
{
  unsigned nsub = offsets[n];
  void* out = malloc(typesize * nsub * width);
  unsigned stride = typesize * width;
  char const* ip = a;
  char* op = out;
  for (unsigned i = 0; i < n; ++i, ip += stride) {
    if (offsets[i] == offsets[i + 1])
      continue;
    memcpy(op, ip, stride);
    op += stride;
  }
  return out;
}

unsigned* ints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* offsets)
{
  return general_subset(sizeof(unsigned), n, width, a, offsets);
}

double* doubles_subset(
    unsigned n,
    unsigned width,
    double const* a,
    unsigned const* offsets)
{
  return general_subset(sizeof(double), n, width, a, offsets);
}

void concat_verts_of_elems(
    unsigned elem_dim,
    unsigned nelems,
    unsigned ngen_elems,
    unsigned const* verts_of_elems,
    unsigned const* offset_of_same_elems,
    unsigned const* verts_of_gen_elems,
    unsigned* nelems_out,
    unsigned** verts_of_elems_out)
{
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned nsame_elems = offset_of_same_elems[nelems];
  unsigned* verts_of_same_elems = ints_subset(nelems, verts_per_elem,
      verts_of_elems, offset_of_same_elems);
  unsigned* out = concat_ints(verts_per_elem,
      verts_of_same_elems, nsame_elems,
      verts_of_gen_elems, ngen_elems);
  free(verts_of_same_elems);
  *nelems_out = nsame_elems + ngen_elems;
  *verts_of_elems_out = out;
}
