#ifndef DERIVE_FACES_H
#define DERIVE_FACES_H

unsigned* derive_faces(
    unsigned nfaces,
    unsigned const* verts_of_elems,
    unsigned const* elems_of_faces,
    unsigned const* elem_face_of_faces);

#endif
