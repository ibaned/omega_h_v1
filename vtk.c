#include "vtk.h"
#include <stdio.h>   // for fprintf, FILE, fclose, fopen, printf
#include "mesh.h"    // for mesh_count, mesh_dim, mesh_find_label
#include "tables.h"  // for the_down_degrees
#include "loop.h"    // for loop_host_malloc
#include <string.h>  // for strlen
#include <assert.h>  // for assert
#include <stdlib.h>  // for atoi

static unsigned const cell_types[4] = {
  1,
  3,
  5,
  10
};

static char const* const type_names[TAG_TYPES] = {
  [TAG_U32] = "UInt32",
  [TAG_F64] = "Float64"
};

static void write_tag(FILE* file, unsigned nents, struct const_tag* tag)
{
  fprintf(file, "<DataArray type=\"%s\" Name=\"%s\""
             " NumberOfComponents=\"%u\" format=\"ascii\">\n",
             type_names[tag->type], tag->name, tag->ncomps);
  switch (tag->type) {
    case TAG_U32: {
      unsigned const* p = tag->data;
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < tag->ncomps; ++j)
          fprintf(file, " %u", *p++);
        fprintf(file, "\n");
      }
    }
    case TAG_F64: {
      double const* p = tag->data;
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < tag->ncomps; ++j)
          fprintf(file, " %e", *p++);
        fprintf(file, "\n");
      }
    }
  }
  fprintf(file, "</DataArray>\n");
}

static char const* types_header =
"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";

/* this function can be a time hog,
 * no fault of our own really, just printf and friends
 * are fairly slow.
 * if you're so inclined, add binary functionality
 * (the VTK format supports it)
 */

void write_vtk(struct mesh* m, char const* filename)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nverts = mesh_count(m, 0);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  FILE* file = fopen(filename, "w");
  fprintf(file, "<VTKFile type=\"UnstructuredGrid\">\n");
  fprintf(file, "<UnstructuredGrid>\n");
  fprintf(file, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n", nverts, nelems);
  fprintf(file, "<Points>\n");
  struct const_tag* coord_tag = mesh_find_tag(m, 0, "coordinates");
  write_tag(file, nverts, coord_tag);
  fprintf(file, "</Points>\n");
  fprintf(file, "<Cells>\n");
  fprintf(file, "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n");
  unsigned down_degree = the_down_degrees[elem_dim][0];
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* p = verts_of_elems + i * down_degree;
    for (unsigned j = 0; j < down_degree; ++j)
      fprintf(file, " %u", p[j]);
    fprintf(file, "\n");
  }
  fprintf(file, "</DataArray>\n");
  fprintf(file, "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n");
  for (unsigned i = 0; i < nelems; ++i)
    fprintf(file, "%u\n", (i + 1) * down_degree);
  fprintf(file, "</DataArray>\n");
  fprintf(file, "%s\n", types_header);
  unsigned type = cell_types[elem_dim];
  for (unsigned i = 0; i < nelems; ++i)
    fprintf(file, "%u\n", type);
  fprintf(file, "</DataArray>\n");
  fprintf(file, "</Cells>\n");
  fprintf(file, "<PointData>\n");
  for (unsigned i = 0; i < mesh_count_tags(m, 0); ++i) {
    struct const_tag* tag = mesh_get_tag(m, 0, i);
    if (tag != coord_tag)
      write_tag(file, nverts, tag);
  }
  fprintf(file, "</PointData>\n");
  fprintf(file, "<CellData>\n");
  for (unsigned i = 0; i < mesh_count_tags(m, mesh_dim(m)); ++i) {
    struct const_tag* tag = mesh_get_tag(m, mesh_dim(m), i);
    write_tag(file, nelems, tag);
  }
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n");
  fprintf(file, "</UnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");
  fclose(file);
}

static char const* the_step_prefix = 0;
static unsigned the_step = 0;

void start_vtk_steps(char const* prefix)
{
  the_step_prefix = prefix;
  the_step = 0;
}

void write_vtk_step(struct mesh* m)
{
  char fname[64];
  sprintf(fname, "%s_%04u.vtu", the_step_prefix, the_step);
  write_vtk(m, fname);
  ++the_step;
}

static unsigned seek_prefix_next(FILE* f,
    char line[], unsigned line_size, char const prefix[])
{
  unsigned pl = (unsigned) strlen(prefix);
  if (!fgets(line, (int) line_size, f))
    return 0;
  return !strncmp(line, prefix, pl);
}

static void seek_prefix(FILE* f,
    char line[], unsigned line_size, char const prefix[])
{
  unsigned pl = (unsigned) strlen(prefix);
  while (fgets(line, (int) line_size, f))
    if (!strncmp(line, prefix, pl))
      return;
  abort();
}

typedef char line_t[1024];

static void read_attrib(char const elem[], char const name[],
    char val[])
{
  char const* pname = strstr(elem, name);
  assert(pname);
  line_t assign;
  strcpy(assign, pname);
  assert(assign[strlen(name) + 1] == '\"');
  char const* pval = strtok(assign + strlen(name) + 2, "\"");
  assert(pval && strlen(pval));
  strcpy(val, pval);
}

static void read_array_name(char const header[], char name[])
{
  read_attrib(header, "Name", name);
}

static enum tag_type read_array_type(char const header[])
{
  line_t val;
  read_attrib(header, "type", val);
  for (unsigned type = 0; type < TAG_TYPES; ++type)
    if (!strcmp(type_names[type], val))
      return type;
  assert(0);
}

static unsigned read_int_attrib(char const header[], char const attrib[])
{
  line_t val;
  read_attrib(header, attrib, val);
  return (unsigned) atoi(val);
}

static unsigned read_array_ncomps(char header[])
{
  return read_int_attrib(header, "NumberOfComponents");
}

static unsigned* read_ints(FILE* f, unsigned n)
{
  unsigned* out = loop_host_malloc(sizeof(unsigned) * n);
  for (unsigned i = 0; i < n; ++i)
    fscanf(f, "%u", &out[i]);
  return out;
}

static double* read_doubles(FILE* f, unsigned n)
{
  double* out = loop_host_malloc(sizeof(double) * n);
  for (unsigned i = 0; i < n; ++i)
    fscanf(f, "%lf", &out[i]);
  return out;
}

static void read_size(FILE* f, unsigned* nverts, unsigned* nelems)
{
  line_t line;
  seek_prefix(f, line, sizeof(line), "<Piece");
  *nverts = read_int_attrib(line, "NumberOfPoints");
  *nelems = read_int_attrib(line, "NumberOfCells");
}

static unsigned read_dimension(FILE* f, unsigned nelems)
{
  assert(nelems);
  line_t line;
  seek_prefix(f, line, sizeof(line), types_header);
  unsigned* types = read_ints(f, nelems);
  unsigned dim;
  for (dim = 0; dim < 4; ++dim)
    if (types[0] == cell_types[dim])
      break;
  assert(dim < 4);
  for (unsigned i = 1; i < nelems; ++i)
    assert(types[i] == cell_types[dim]);
  loop_host_free(types);
  return dim;
}

static unsigned read_mesh_array(FILE* f, struct mesh* m,
    unsigned dim)
{
  line_t line;
  if (!seek_prefix_next(f, line, sizeof(line), "<DataArray"))
    return 0;
  enum tag_type type = read_array_type(line);
  line_t name;
  read_array_name(line, name);
  unsigned ncomps = read_array_ncomps(line);
  void* data;
  switch (type) {
    case TAG_U32: data = read_ints(f, mesh_count(m, dim) * ncomps);
    case TAG_F64: data = read_doubles(f, mesh_count(m, dim) * ncomps);
  }
  mesh_add_tag(m, dim, type, name, ncomps, data);
  seek_prefix(f, line, sizeof(line), "</DataArray");
  return 1;
}

static void read_verts(FILE* f, struct mesh* m)
{
  line_t line;
  seek_prefix(f, line, sizeof(line), "<Points");
  unsigned ok = read_mesh_array(f, m, 0);
  assert(ok);
}

static void read_elems(FILE* f, struct mesh* m, unsigned nelems)
{
  line_t line;
  seek_prefix(f, line, sizeof(line), "<Cells");
  seek_prefix(f, line, sizeof(line), "<DataArray");
  line_t name;
  read_array_name(line, name);
  assert(!strcmp("connectivity", name));
  unsigned dim = mesh_dim(m);
  unsigned verts_per_elem = the_down_degrees[dim][0];
  unsigned* data = read_ints(f, nelems * verts_per_elem);
  mesh_set_ents(m, dim, nelems, data);
}

static struct mesh* read_vtk_mesh(FILE* f)
{
  unsigned nverts, nelems;
  read_size(f, &nverts, &nelems);
  if (!nelems) {
    fclose(f);
    return new_mesh(0);
  }
  unsigned dim = read_dimension(f, nelems);
  struct mesh* m = new_mesh(dim);
  mesh_set_ents(m, 0, nverts, 0);
  rewind(f);
  read_verts(f, m);
  read_elems(f, m, nelems);
  return m;
}

static void read_vtk_fields(FILE* f, struct mesh* m)
{
  line_t line;
  seek_prefix(f, line, sizeof(line), "<PointData");
  while(read_mesh_array(f, m, 0));
  seek_prefix(f, line, sizeof(line), "<CellData");
  while(read_mesh_array(f, m, mesh_dim(m)));
}

struct mesh* read_vtk(char const* filename)
{
  FILE* f = fopen(filename, "r");
  assert(f);
  struct mesh* m = read_vtk_mesh(f);
  read_vtk_fields(f, m);
  fclose(f);
  return m;
}
