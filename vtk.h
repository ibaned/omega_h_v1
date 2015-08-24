#ifndef VTK_H
#define VTK_H

#include "mesh.h"

void write_vtk(struct mesh* m, char const* filename);

void start_vtk_steps(char const* prefix);
void write_vtk_step(struct mesh* m);

#endif
