#include <stdio.h>

#include "arrays.h"
#include "comm.h"
#include "loop.h"
#include "mesh.h"

int main()
{
  comm_init();
  /* make a 2-triangle mesh with only
     the triangle to vertex graph */
  struct mesh* m = new_box_mesh(2);
  /* derive the edges, resulting in at least
     the edge to vertex graph */
  mesh_ask_down(m, 1, 0);
  /* use reflect_down to derive the
     triangle to edge graph */
  mesh_ask_down(m, 2, 1);
  unsigned* host_graph = uints_to_host(
      mesh_ask_down(m, 2, 1),
      mesh_count(m, 2) * 3);
  for (unsigned i = 0; i < mesh_count(m, 2); ++i) {
    for (unsigned j = 0; j < 3; ++j)
      printf("%u ", host_graph[i * 3 + j]);
    printf("\n");
  }
  loop_free(host_graph);
  free_mesh(m);
  comm_fini();
}
