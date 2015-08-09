#ifndef CONCAT_H
#define CONCAT_H

unsigned* ints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* offsets);

double* doubles_subset(
    unsigned n,
    unsigned width,
    double const* a,
    unsigned const* offsets);

unsigned* concat_ints(
    unsigned narrs,
    unsigned width,
    unsigned const* sizes,
    unsigned const* const* arrs);

double* concat_doubles(
    unsigned narrs,
    unsigned width,
    unsigned const* sizes,
    double const* const* arrs);

void concat_verts_of_elems(
    unsigned elem_dim,
    unsigned nelems,
    unsigned ngen_elems,
    unsigned const* verts_of_elems,
    unsigned const* offset_of_same_elems,
    unsigned const* verts_of_gen_elems,
    unsigned* nelems_out,
    unsigned** verts_of_elems_out);

#endif
