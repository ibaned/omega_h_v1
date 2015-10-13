#include "derive_faces.h"

#include "loop.h"
#include "tables.h"

unsigned* derive_faces(
    unsigned nfaces,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_faces,
    unsigned const* elem_face_of_faces)
{
  unsigned* verts_of_faces = LOOP_MALLOC(unsigned, nfaces * 3);
  unsigned const* const* elem_verts_of_faces =
    the_canonical_orders[3][2][0];
  for (unsigned i = 0; i < nfaces; ++i) {
    unsigned elem = elems_of_faces[i * 2];
    unsigned const* verts_of_elem = verts_of_elems + elem * 4;
    unsigned elem_face = elem_face_of_faces[i];
    unsigned const* elem_verts_of_face = elem_verts_of_faces[elem_face];
    unsigned* verts_of_face = verts_of_faces + i * 3;
    for (unsigned j = 0; j < 3; ++j) {
      unsigned elem_vert_of_face = elem_verts_of_face[j];
      unsigned vert = verts_of_elem[elem_vert_of_face];
      verts_of_face[j] = vert;
    }
  }
  return verts_of_faces;
}
