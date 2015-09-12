#include "tables.h"
#include "up_from_down.h"
#include "reflect_down.h"
#include "bridge_graph.h"
#include "derive_faces.h"
#include "loop.h"
#include <stdio.h>

int main()
{
  unsigned elem_dim = 3;
  unsigned nelems;
  unsigned nverts;
  unsigned* verts_of_elems;
  get_box_copy(elem_dim, &nelems, &nverts, &verts_of_elems, 0);
  unsigned* elems_of_verts_offsets;
  unsigned* elems_of_verts;
  up_from_down(elem_dim, 0, nelems, nverts, verts_of_elems,
      &elems_of_verts_offsets, &elems_of_verts, 0);
  unsigned* elems_of_elems = get_dual(elem_dim, nelems,
      verts_of_elems, elems_of_verts_offsets, elems_of_verts);
  loop_free(elems_of_verts_offsets);
  loop_free(elems_of_verts);
  unsigned nfaces = 0;
  unsigned* elems_of_faces = 0;
  unsigned* elem_face_of_faces = 0;
  bridge_dual_graph(elem_dim, nelems, elems_of_elems,
      &nfaces, &elems_of_faces, &elem_face_of_faces);
  loop_free(elems_of_elems);
  for (unsigned i = 0; i < nfaces; ++i) {
    printf("face %u: ", i);
    printf("elems %u %u, ", elems_of_faces[i * 2 + 0],
                            elems_of_faces[i * 2 + 1]);
    printf("face %u of elem %u\n",
        elem_face_of_faces[i],
        elems_of_faces[i * 2]);
  }
  unsigned* verts_of_faces = derive_faces(nfaces, verts_of_elems,
      elems_of_faces, elem_face_of_faces);
  loop_free(verts_of_elems);
  loop_free(elems_of_faces);
  loop_free(elem_face_of_faces);
  printf("\n");
  for (unsigned i = 0; i < nfaces; ++i) {
    printf("face %u:", i);
    for (unsigned j = 0; j < 3; ++j)
      printf(" %u", verts_of_faces[i * 3 + j]);
    printf("\n");
  }
  loop_free(verts_of_faces);
}
