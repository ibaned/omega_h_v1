#ifndef GLOBAL_H
#define GLOBAL_H

struct mesh;

void mesh_number_simply(struct mesh* m);
unsigned long* globalize_offsets(unsigned* local, unsigned n);
void global_to_linpart(unsigned long* global, unsigned n,
    unsigned long total, unsigned nparts,
    unsigned** p_part, unsigned** p_local);
void categorize_by_part(unsigned const* parts, unsigned n,
    unsigned** p_cats, unsigned** p_cat_indices,
    unsigned* p_ncats,
    unsigned** p_cat_parts, unsigned** p_cat_counts);
unsigned* sort_uints_by_category(
    unsigned const* a,
    unsigned n,
    unsigned const* cats,
    unsigned const* indices,
    unsigned const* cat_offsets);

#endif
