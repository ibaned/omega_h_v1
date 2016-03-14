#ifndef GMSH_IO_HPP
#define GMSH_IO_HPP

namespace omega_h {

struct mesh;

struct mesh* read_msh(char const* filename);

}

#endif
