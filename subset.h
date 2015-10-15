#ifndef SUBSET_H
#define SUBSET_H

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

void vert_tags_subset(struct mesh* in, struct mesh* out,
    unsigned const* offsets);

struct mesh* subset_mesh(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* offsets);

#endif
