#include "concat.h"
#include "tables.h"
#include "ints.h"
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
  unsigned b_bytes = nb * width * typesize;
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

void graph_subset(
    unsigned n,
    unsigned const* set_offsets,
    unsigned const* graph_offsets,
    unsigned const* graph_adj,
    unsigned** graph_offsets_out,
    unsigned** graph_adj_out)
{
  unsigned n_out = set_offsets[n];
  unsigned* degrees = ints_unscan(set_offsets, n);
  unsigned* degrees_out = ints_subset(n, 1, degrees, set_offsets);
  free(degrees);
  unsigned* new_offsets = ints_exscan(degrees_out, n_out);
  free(degrees_out);
  unsigned nadj_out = new_offsets[n_out];
  unsigned* new_adj = malloc(sizeof(unsigned) * nadj_out);
  for (unsigned i = 0; i < n; ++i) {
    if (set_offsets[i] == set_offsets[i + 1])
      continue;
    unsigned new_i = set_offsets[i];
    unsigned first_adj = graph_offsets[i];
    unsigned end_adj = graph_offsets[i];
    unsigned first_new_adj = new_offsets[new_i];
    for (unsigned j = first_adj; j < end_adj; ++j)
      new_adj[(j - first_adj) + first_new_adj] = graph_adj[j];
  }
  *graph_offsets_out = new_offsets;
  *graph_adj_out = new_adj;
}
