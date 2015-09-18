#ifndef VTK_H
#define VTK_H

struct mesh;

void write_vtk(struct mesh* m, char const* filename);

void start_vtk_steps(char const* prefix);
void write_vtk_step(struct mesh* m);

struct mesh* read_vtk(char const* filename);

#endif
