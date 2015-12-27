#ifndef SUBSET_H
#define SUBSET_H

void uchars_subset_into(
    unsigned n,
    unsigned width,
    unsigned char const* a,
    unsigned const* offsets,
    unsigned char* out);

void uints_subset_into(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* offsets,
    unsigned* out);

void ulongs_subset_into(
    unsigned n,
    unsigned width,
    unsigned long const* a,
    unsigned const* offsets,
    unsigned long* out);

void doubles_subset_into(
    unsigned n,
    unsigned width,
    double const* a,
    unsigned const* offsets,
    double* into);

unsigned char* uchars_subset(
    unsigned n,
    unsigned width,
    unsigned char const* a,
    unsigned const* offsets);

unsigned* uints_subset(
    unsigned n,
    unsigned width,
    unsigned const* a,
    unsigned const* offsets);

unsigned long* ulongs_subset(
    unsigned n,
    unsigned width,
    unsigned long const* a,
    unsigned const* offsets);

double* doubles_subset(
    unsigned n,
    unsigned width,
    double const* a,
    unsigned const* offsets);

struct mesh;

void tags_subset(struct mesh* in, struct mesh* out,
    unsigned dim, unsigned const* offsets);

struct mesh* subset_mesh(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* offsets);

void subset_verts_of_doms(
    struct mesh* m,
    unsigned dom_dim,
    unsigned const* offset_of_doms,
    unsigned* verts_of_prods);

#endif
