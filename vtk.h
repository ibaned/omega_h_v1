#ifndef VTK_H
#define VTK_H

struct mesh;
struct cloud;

enum vtk_format {
  VTK_ASCII,
  VTK_BINARY
};

#define VTK_FORMATS 2

void write_vtu(struct mesh* m, char const* filename);

void start_vtk_steps(char const* prefix);
void write_vtk_step(struct mesh* m);

struct mesh* read_vtu(char const* filename);

void write_vtu_cloud(struct cloud* c, char const* filename);
struct cloud* read_vtu_cloud(char const* filename);

void write_pvtu(struct mesh* m, char const* filename,
    unsigned npieces, unsigned nghost_levels);

struct mesh* read_parallel_vtu(char const* inpath);
void write_parallel_vtu(struct mesh* m, char const* outpath);

#endif
