#include "node_ele_io.h"
#include "mesh.h"
#include "loop.h"
#include "field.h"
#include "label.h"
#include <stdio.h>
#include <assert.h>

struct mesh* read_dot_node(char const* filename)
{
  FILE* f = fopen(filename, "r");
  unsigned nnodes, dim, nattrib, nbdrys;
  fscanf(f, "%u %u %u %u", &nnodes, &dim, &nattrib, &nbdrys);
  struct mesh* m = new_mesh(dim);
  double* coords = loop_host_malloc(sizeof(double) * nnodes * 3);
  assert(nbdrys == 0 || nbdrys == 1);
  double* attrib = 0;
  if (nattrib)
    attrib = loop_host_malloc(sizeof(double) * nnodes * nattrib);
  unsigned* bdry = 0;
  if (nbdrys)
    bdry = loop_host_malloc(sizeof(unsigned) * nnodes);
  for (unsigned i = 0; i < nnodes; ++i) {
    fscanf(f, "%*u");
    unsigned j;
    for (j = 0; j < dim; ++j)
      fscanf(f, "%lf", coords + i * 3 + j);
    for (; j < 3; ++j)
      coords[i * 3 + j] = 0;
    for (j = 0; j < nattrib; ++j)
      fscanf(f, "%lf", attrib + i * nattrib + j);
    if (nbdrys)
      fscanf(f, "%u", bdry + i);
  }
  fclose(f);
  mesh_set_ents(m, 0, nnodes, 0);
  mesh_add_nodal_field(m, "coordinates", 3, coords);
  if (nattrib)
    mesh_add_nodal_field(m, "attributes", nattrib, attrib);
  if (nbdrys)
    mesh_add_nodal_label(m, "boundary", bdry);
  return m;
}

void write_dot_node(struct mesh* m, char const* filename)
{
  double const* coords =
    mesh_find_nodal_field(m, "coordinates")->data;
  struct const_field* attribf =
    mesh_find_nodal_field(m, "attributes");
  struct const_label* bdryl =
    mesh_find_nodal_label(m, "boundary");
  unsigned nnodes = mesh_count(m, 0);
  unsigned dim = mesh_dim(m);
  unsigned nattrib = 0;
  if (attribf)
    nattrib = attribf->ncomps;
  unsigned nbdrys = 0;
  if (bdryl)
    nbdrys = 1;
  FILE* f = fopen(filename, "w");
  fprintf(f, "%u %u %u %u\n", nnodes, dim, nattrib, nbdrys);
  for (unsigned i = 0; i < nnodes; ++i) {
    fprintf(f, " %u", i + 1);
    for (unsigned j = 0; j < dim; ++j)
      fprintf(f, " %f", coords[i * 3 + j]);
    for (unsigned j = 0; j < nattrib; ++j)
      fprintf(f, " %f", attribf->data[i * nattrib + j]);
    if (nbdrys)
      fprintf(f, " %u", bdryl->data[i]);
    fprintf(f, "\n");
  }
  fclose(f);
}

void read_dot_ele(struct mesh* m, char const* filename)
{
  unsigned nelems, verts_per_elem, nattrib;
  FILE* f = fopen(filename, "r");
  fscanf(f, "%u %u %u", &nelems, &verts_per_elem, &nattrib);
  assert(verts_per_elem >= 2);
  unsigned elem_dim = verts_per_elem - 1;
  assert(elem_dim <= mesh_dim(m));
  unsigned* verts_of_elems = loop_host_malloc(
      sizeof(unsigned) * nelems * verts_per_elem);
  double* attrib = 0;
  if (nattrib)
    attrib = loop_host_malloc(sizeof(double) * nelems * nattrib);
  for (unsigned i = 0; i < nelems; ++i) {
    fscanf(f, "%*u");
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned v;
      fscanf(f, "%u", &v);
      verts_of_elems[i * verts_per_elem + j] = v - 1;
    }
    for (unsigned j = 0; j < nattrib; ++j)
      fscanf(f, "%lf", attrib + i * nattrib + j);
  }
  fclose(f);
  mesh_set_ents(m, elem_dim, nelems, verts_of_elems);
  if (nattrib)
    mesh_add_elem_field(m, "attributes", nattrib, attrib);
}

void write_dot_ele(struct mesh* m, char const* filename)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned verts_per_elem = elem_dim + 1;
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned nattrib = 0;
  struct const_field* attribf = mesh_find_elem_field(m, "attributes");
  if (attribf)
    nattrib = attribf->ncomps;
  FILE* f = fopen(filename, "w");
  fprintf(f, "%u %u %u\n", nelems, verts_per_elem, nattrib);
  for (unsigned i = 0; i < nelems; ++i) {
    fprintf(f, " %u", i + 1);
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned v = verts_of_elems[i * verts_per_elem + j];
      fprintf(f, " %u", v + 1);
    }
    for (unsigned j = 0; j < nattrib; ++j)
      fprintf(f, " %f", attribf->data[i * nattrib + j]);
    fprintf(f, "\n");
  }
  fclose(f);
}
