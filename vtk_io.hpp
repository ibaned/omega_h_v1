#ifndef VTK_IO_HPP
#define VTK_IO_HPP

namespace omega_h {

struct mesh;

enum vtk_format {
  VTK_ASCII,
  VTK_BINARY
};

#define VTK_FORMATS 2

struct mesh* read_mesh_vtk(char const* filename);
void write_mesh_vtk(struct mesh* m, char const* filename);
void write_mesh_vtk_opts(struct mesh* m, char const* outpath,
    enum vtk_format fmt);

void start_vtk_steps(char const* prefix);
void write_vtk_step(struct mesh* m);

}

#endif
