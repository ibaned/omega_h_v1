#ifndef INDSET_H
#define INDSET_H

/* given an undirected graph (usually obtained from get_star),
 * an initial filtered subset of the vertices,
 * and a scalar value (goodness) at the vertices,
 * returns a maximal independent set of the subgraph
 * induced by the filter, whose members are
 * preferrably those with locally high goodness values.
 * In order for this function to run efficiently,
 * there should not exist long paths with monotonically
 * decreasing goodness.
 */

unsigned* find_indset(
    unsigned nverts,
    unsigned const* offsets,
    unsigned const* adj,
    unsigned const* filter,
    double const* goodness);

struct mesh;

unsigned* mesh_find_indset(struct mesh* m, unsigned ent_dim,
    unsigned const* candidates, double const* qualities);

unsigned* mesh_indset_offsets(struct mesh* m, unsigned ent_dim,
    unsigned const* candidates, double const* qualities);

#endif
