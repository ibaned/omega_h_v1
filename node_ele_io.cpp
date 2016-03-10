#include "node_ele_io.h"

#include <assert.h>
#include <stdio.h>

#include "files.h"
#include "loop.h"
#include "mesh.h"
#include "tag.h"

struct mesh* read_dot_node(char const* filename)
{
  FILE* f = safe_fopen(filename, "r");
  unsigned nnodes, dim, nattrib, nbdrys;
  safe_scanf(f, 4, "%u %u %u %u", &nnodes, &dim, &nattrib, &nbdrys);
  assert(nnodes < 100 * 1000 * 1000);
  assert(nattrib < 100);
  struct mesh* m = new_mesh(dim, MESH_REDUCED, 0);
  double* coords = LOOP_HOST_MALLOC(double, nnodes * 3);
  assert(nbdrys == 0 || nbdrys == 1);
  double* attrib = 0;
  if (nattrib)
    attrib = LOOP_HOST_MALLOC(double, nnodes * nattrib);
  unsigned* bdry = 0;
  if (nbdrys)
    bdry = LOOP_HOST_MALLOC(unsigned, nnodes);
  for (unsigned i = 0; i < nnodes; ++i) {
    safe_scanf(f, 0, "%*u");
    unsigned j;
    for (j = 0; j < dim; ++j)
      safe_scanf(f, 1, "%lf", coords + i * 3 + j);
    for (; j < 3; ++j)
      coords[i * 3 + j] = 0;
    for (j = 0; j < nattrib; ++j)
      safe_scanf(f, 1, "%lf", attrib + i * nattrib + j);
    if (nbdrys)
      safe_scanf(f, 1, "%u", bdry + i);
  }
  fclose(f);
  mesh_set_ents(m, 0, nnodes, 0);
  mesh_add_tag(m, 0, TAG_F64, "coordinates", 3, coords);
  if (nattrib)
    mesh_add_tag(m, 0, TAG_F64, "attributes", nattrib, attrib);
  if (nbdrys)
    mesh_add_tag(m, 0, TAG_U32, "boundary", 1, bdry);
  return m;
}

void write_dot_node(struct mesh* m, char const* filename)
{
  double const* coords =
    mesh_find_tag(m, 0, "coordinates")->d.f64;
  struct const_tag* attribf =
    mesh_find_tag(m, 0, "attributes");
  struct const_tag* bdryl =
    mesh_find_tag(m, 0, "boundary");
  unsigned nnodes = mesh_count(m, 0);
  unsigned dim = mesh_dim(m);
  unsigned nattrib = 0;
  double const* attrib = 0;
  if (attribf) {
    nattrib = attribf->ncomps;
    attrib = attribf->d.f64;
  }
  unsigned nbdrys = 0;
  unsigned const* bdry = 0;
  if (bdryl) {
    nbdrys = 1;
    bdry = bdryl->d.u32;
  }
  FILE* f = safe_fopen(filename, "w");
  fprintf(f, "%u %u %u %u\n", nnodes, dim, nattrib, nbdrys);
  for (unsigned i = 0; i < nnodes; ++i) {
    fprintf(f, " %u", i + 1);
    for (unsigned j = 0; j < dim; ++j)
      fprintf(f, " %f", coords[i * 3 + j]);
    for (unsigned j = 0; j < nattrib; ++j)
      fprintf(f, " %f", attrib[i * nattrib + j]);
    if (nbdrys)
      fprintf(f, " %u", bdry[i]);
    fprintf(f, "\n");
  }
  fclose(f);
}

void read_dot_ele(struct mesh* m, char const* filename)
{
  unsigned nelems, verts_per_elem, nattrib;
  FILE* f = safe_fopen(filename, "r");
  safe_scanf(f, 3, "%u %u %u", &nelems, &verts_per_elem, &nattrib);
  assert(nelems < 100 * 1000 * 1000);
  assert(nattrib < 100);
  assert(verts_per_elem >= 2);
  assert(verts_per_elem <= 4);
  unsigned elem_dim = verts_per_elem - 1;
  assert(elem_dim <= mesh_dim(m));
  unsigned* verts_of_elems = LOOP_HOST_MALLOC(unsigned,
      nelems * verts_per_elem);
  double* attrib = 0;
  if (nattrib)
    attrib = LOOP_HOST_MALLOC(double, nelems * nattrib);
  for (unsigned i = 0; i < nelems; ++i) {
    safe_scanf(f, 0, "%*u");
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned v;
      safe_scanf(f, 1, "%u", &v);
      verts_of_elems[i * verts_per_elem + j] = v - 1;
    }
    for (unsigned j = 0; j < nattrib; ++j)
      safe_scanf(f, 1, "%lf", attrib + i * nattrib + j);
  }
  fclose(f);
  mesh_set_ents(m, elem_dim, nelems, verts_of_elems);
  if (nattrib)
    mesh_add_tag(m, elem_dim, TAG_F64, "attributes", nattrib, attrib);
}

void write_dot_ele(struct mesh* m, char const* filename)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned verts_per_elem = elem_dim + 1;
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned nattrib = 0;
  struct const_tag* attribf = mesh_find_tag(m, elem_dim, "attributes");
  double const* attrib = 0;
  if (attribf) {
    nattrib = attribf->ncomps;
    attrib = attribf->d.f64;
  }
  FILE* f = safe_fopen(filename, "w");
  fprintf(f, "%u %u %u\n", nelems, verts_per_elem, nattrib);
  for (unsigned i = 0; i < nelems; ++i) {
    fprintf(f, " %u", i + 1);
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned v = verts_of_elems[i * verts_per_elem + j];
      fprintf(f, " %u", v + 1);
    }
    for (unsigned j = 0; j < nattrib; ++j)
      fprintf(f, " %f", attrib[i * nattrib + j]);
    fprintf(f, "\n");
  }
  fclose(f);
}
