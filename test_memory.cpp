#include <assert.h>
#include <stdio.h>

#include "loop.hpp"
#include "vtk_io.hpp"
#include "mesh.hpp"

int main(int argc, char** argv)
{
  assert(argc == 2);
  printf("baseline %lu\n", loop_host_memory());
  struct mesh* m = read_mesh_vtk(argv[1]);
  printf("regions to vertices %lu\n", loop_host_memory());
  printf("regions to vertices %f (per element)\n",
      ((double)(loop_host_memory())/((double)(mesh_count(m, mesh_dim(m))))));
  mesh_ask_up(m, 0, mesh_dim(m));
  printf("regions <-> vertices %lu\n", loop_host_memory());
  printf("regions <-> vertices %f (per element)\n",
      ((double)(loop_host_memory())/((double)(mesh_count(m, mesh_dim(m))))));
  mesh_count(m, 1);
  printf("with edges %lu\n", loop_host_memory());
  printf("with edges %f (per element)\n",
      ((double)(loop_host_memory())/((double)(mesh_count(m, mesh_dim(m))))));
  mesh_count(m, 2);
  printf("with edges and faces %lu\n", loop_host_memory());
  printf("with edges and faces %f (per element)\n",
      ((double)(loop_host_memory())/((double)(mesh_count(m, mesh_dim(m))))));
  free_mesh(m);
  printf("after freeing %lu\n", loop_host_memory());
}
