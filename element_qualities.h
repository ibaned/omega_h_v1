#ifndef ELEMENT_QUALITIES_H
#define ELEMENT_QUALITIES_H

double* element_qualities(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords);

double min_element_quality(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords);

#endif
