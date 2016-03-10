#include <assert.h>
#include <stdio.h>

#include "comm.hpp"
#include "mesh.hpp"
#include "vtk_io.hpp"

static void write_graph(char const* filename,
    unsigned n,
    unsigned const* off,
    unsigned const* adj,
    double const* xyz)
{
  FILE* f = fopen(filename, "w");
  fprintf(f, "%u\n", n);
  for (unsigned i = 0; i < n; ++i)
    fprintf(f, "%u\n", off[i + 1]);
  for (unsigned i = 0; i < off[n]; ++i)
    fprintf(f, "%u\n", adj[i]);
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < 3; ++j)
      fprintf(f, " %f", xyz[i * 3 + j]);
    fprintf(f, "\n");
  }
  fclose(f);
}

int main(int argc, char** argv)
{
  assert(argc == 3);
  comm_init();
  struct mesh* m = read_mesh_vtk(argv[1]);
  write_graph(argv[2],
      mesh_count(m, 0),
      mesh_ask_star(m, 0, 1)->offsets,
      mesh_ask_star(m, 0, 1)->adj,
      mesh_find_tag(m, 0, "coordinates")->d.f64);
  free_mesh(m);
  comm_fini();
}

