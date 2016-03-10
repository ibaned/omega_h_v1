#ifndef REFINE_TOPOLOGY_H
#define REFINE_TOPOLOGY_H

/* this is the central function involved in
 * mesh refinement.
 * there are three variable dimensions:
 *   the "domain" dimension (interiors are being filled)
 *   the "source" dimension (new vertices at the center)
 *   the "product" dimension (filler for domain interior)
 * this function generates the "product" entities
 * which fill the open domain of the "domain" entities
 * as a result of certain "source" entities being
 * split in the middle by a new vertex.
 *
 * both domains and products are represented in terms
 * of their vertices
 */

void refine_topology(
    unsigned dom_dim,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned ndoms,
    unsigned const* verts_of_doms,
    unsigned const* offset_of_doms,
    unsigned const* direction_of_doms,
    unsigned const* vert_of_doms,
    unsigned* verts_of_prods);

unsigned get_prods_per_dom(
    unsigned dom_dim,
    unsigned src_dim,
    unsigned prod_dim);

struct mesh;

void mesh_refine_topology(struct mesh* m,
    unsigned dom_dim,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned const* offset_of_doms,
    unsigned const* vert_of_doms,
    unsigned const* direction_of_doms,
    unsigned* verts_of_prods);

#endif
