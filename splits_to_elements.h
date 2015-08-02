#ifndef SPLITS_TO_ELEMENTS_H
#define SPLITS_TO_ELEMENTS_H

struct splits_to_elements {
  unsigned* elem_split_offset;
  unsigned* elem_split_vert;
  unsigned* elem_split_direction;
};

struct splits_to_elements project_splits_to_elements(
    unsigned elem_dim,
    unsigned split_dim,
    unsigned nelem,
    unsigned const* elem_splits,
    unsigned const* split_offsets,
    unsigned const* split_new_verts);

#endif
