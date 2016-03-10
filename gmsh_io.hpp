#ifndef GMSH_IO_H
#define GMSH_IO_H

struct mesh;

struct mesh* read_msh(char const* filename);

#endif
