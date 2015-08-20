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
    unsigned width,
    unsigned const a[],
    unsigned na,
    unsigned const b[],
    unsigned nb);

double* concat_doubles(
    unsigned width,
    double const a[],
    unsigned na,
    double const b[],
    unsigned nb);

void concat_verts_of_elems(
    unsigned elem_dim,
    unsigned nelems,
    unsigned ngen_elems,
    unsigned const* verts_of_elems,
    unsigned const* offset_of_same_elems,
    unsigned const* verts_of_gen_elems,
    unsigned* nelems_out,
    unsigned** verts_of_elems_out);

void graph_subset(
    unsigned n,
    unsigned const* set_offsets,
    unsigned const* graph_offsets,
    unsigned const* graph_adj,
    unsigned** graph_offsets_out,
    unsigned** graph_adj_out);

void concat_graph(
    unsigned na,
    unsigned const* a_offsets,
    unsigned const* a_adj,
    unsigned nb,
    unsigned const* b_offsets,
    unsigned const* b_adj,
    unsigned** p_offsets,
    unsigned** p_adj);

#endif
