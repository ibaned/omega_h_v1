#include "rv_mesh.h"
#include "tables.h"
#include <stdlib.h>
#include <string.h>

struct rv_mesh new_box_rv_mesh(unsigned dim)
{
  struct rv_mesh m;
  m.elem_dim = dim;
  m.nelems = the_box_nelems[dim];
  m.nverts = the_box_nverts[dim];
  unsigned verts_per_elem = the_down_degrees[dim][0];
  unsigned bytes = sizeof(unsigned) * m.nelems * verts_per_elem;
  m.verts_of_elems = malloc(bytes);
  memcpy(m.verts_of_elems, the_box_conns[dim], bytes);
  bytes = sizeof(double) * m.nverts * 3;
  m.coords = malloc(bytes);
  memcpy(m.coords, the_box_coords[dim], bytes);
  return m;
}

void free_rv_mesh(struct rv_mesh m)
{
  free(m.verts_of_elems);
  free(m.coords);
}
