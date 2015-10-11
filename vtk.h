#ifndef VTK_H
#define VTK_H

struct mesh;
struct cloud;

void write_vtk(struct mesh* m, char const* filename);

void start_vtk_steps(char const* prefix);
void write_vtk_step(struct mesh* m);

struct mesh* read_vtk(char const* filename);

void write_vtk_cloud(struct cloud* c, char const* filename);
struct cloud* read_vtk_cloud(char const* filename);

void write_pvtu(struct mesh* m, char const* filename,
    unsigned npieces, unsigned nghost_levels);

#endif
