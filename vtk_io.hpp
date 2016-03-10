#ifndef VTK_IO_HPP
#define VTK_IO_HPP

struct mesh;
struct cloud;

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

void write_vtu_cloud_opts(struct cloud* c, char const* filename,
    enum vtk_format fmt);
void write_vtu_cloud(struct cloud* c, char const* filename);
struct cloud* read_vtu_cloud(char const* filename);

void write_pvtu_cloud(struct cloud* c, char const* filename,
    unsigned npieces);

struct cloud* read_parallel_vtu_cloud(char const* inpath);
void write_parallel_vtu_cloud(struct cloud* c, char const* outpath);

#endif
