#include "vtk.h"
#include "tables.h"
#include <stdio.h>

static unsigned const vtk_types[4] = {
  1,
  3,
  5,
  10
};

void write_vtk(
    char const* filename,
    unsigned elem_dim,
    unsigned nelem,
    unsigned nvert,
    unsigned const* elem_to_vert,
    double const* coords)
{
  FILE* f = fopen(filename, "w");
  fprintf(f, "<VTKFile type=\"UnstructuredGrid\" version=\"3.0\">\n");
  fprintf(f, "<UnstructuredGrid>\n");
  fprintf(f, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\">\n", nvert, nelem);
  fprintf(f, "<Points>\n");
  fprintf(f, "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (unsigned i = 0; i < nvert; ++i) {
    double const* p = coords + i * 3;
    fprintf(f, "%f %f %f\n", p[0], p[1], p[2]);
  }
  fprintf(f, "</DataArray>\n");
  fprintf(f, "</Points>\n");
  fprintf(f, "<Cells>\n");
  fprintf(f, "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n");
  unsigned down_degree = the_down_degrees[elem_dim][0];
  for (unsigned i = 0; i < nelem; ++i) {
    unsigned const* p = elem_to_vert + i * down_degree;
    for (unsigned j = 0; j < down_degree; ++j)
      fprintf(f, " %u", p[j]);
    fprintf(f, "\n");
  }
  fprintf(f, "</DataArray>\n");
  fprintf(f, "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n");
  for (unsigned i = 0; i < nelem; ++i)
    fprintf(f, "%u\n", (i + 1) * down_degree);
  fprintf(f, "</DataArray>\n");
  fprintf(f, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  unsigned type = vtk_types[elem_dim];
  for (unsigned i = 0; i < nelem; ++i)
    fprintf(f, "%u\n", type);
  fprintf(f, "</DataArray>\n");
  fprintf(f, "</Cells>\n");
  fprintf(f, "</Piece>");
  fprintf(f, "</UnstructuredGrid>\n");
  fprintf(f, "</VTKFile>\n");
  fclose(f);
}
