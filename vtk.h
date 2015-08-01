#ifndef VTK_H
#define VTK_H

void write_vtk(
   char const* filename,
   unsigned elem_dim,
   unsigned nelem,
   unsigned nvert,
   unsigned const* elem_to_vert,
   double const* coords);

#endif
