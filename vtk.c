#include "vtk.h"
#include "tables.h"
#include <stdio.h>

static unsigned const vtk_types[4] = {
  1,
  3,
  5,
  10
};

/* this function can be a time hog,
 * no fault of our own really, just printf and friends
 * are fairly slow.
 * if you're so inclined, add binary functionality
 * (the VTK format supports it)
 */

static void write_nodal_field(FILE* file, struct mesh* m, struct const_field* field)
{
  fprintf(file, "<DataArray type=\"Float64\" Name=\"%s\""
             " NumberOfComponents=\"%u\" format=\"ascii\">\n",
             field->name, field->ncomps);
  unsigned nverts = mesh_count(m, 0);
  double const* p = field->data;
  for (unsigned i = 0; i < nverts; ++i) {
    for (unsigned j = 0; j < field->ncomps; ++j)
      fprintf(file, " %f", *p++);
    fprintf(file, "\n");
  }
  fprintf(file, "</DataArray>\n");
}

static void write_elem_field(FILE* file, struct mesh* m, struct const_field* field)
{
  fprintf(file, "<DataArray type=\"Float64\" Name=\"%s\""
             " NumberOfComponents=\"%u\" format=\"ascii\">\n",
             field->name, field->ncomps);
  unsigned nverts = mesh_count(m, mesh_dim(m));
  double const* p = field->data;
  for (unsigned i = 0; i < nverts; ++i) {
    for (unsigned j = 0; j < field->ncomps; ++j)
      fprintf(file, " %f", *p++);
    fprintf(file, "\n");
  }
  fprintf(file, "</DataArray>\n");
}

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
  struct const_field* coord_field = mesh_find_nodal_field(m, "coordinates");
  write_nodal_field(file, m, coord_field);
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
  fprintf(file, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  unsigned type = vtk_types[elem_dim];
  for (unsigned i = 0; i < nelems; ++i)
    fprintf(file, "%u\n", type);
  fprintf(file, "</DataArray>\n");
  fprintf(file, "</Cells>\n");
  fprintf(file, "<PointData>\n");
  if (mesh_find_nodal_label(m, "class_dim")) {
    unsigned const* class_dim = mesh_find_nodal_label(m, "class_dim")->data;
    fprintf(file, "<DataArray type=\"UInt32\" Name=\"class_dim\""
               " NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (unsigned i = 0; i < nverts; ++i)
      fprintf(file, "%u\n", class_dim[i]);
    fprintf(file, "</DataArray>\n");
  }
  for (unsigned i = 0; i < mesh_count_nodal_fields(m); ++i) {
    struct const_field* field = mesh_get_nodal_field(m, i);
    if (field != coord_field)
      write_nodal_field(file, m, field);
  }
  fprintf(file, "</PointData>\n");
  fprintf(file, "<CellData>\n");
  for (unsigned i = 0; i < mesh_count_elem_fields(m); ++i) {
    struct const_field* field = mesh_get_elem_field(m, i);
    if (field != coord_field)
      write_elem_field(file, m, field);
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
