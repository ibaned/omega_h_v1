#ifndef REFINE_TOPOLOGY_H
#define REFINE_TOPOLOGY_H

/* this is the central function involved in
 * mesh refinement.
 * there are three variable dimensions:
 *   the "element" dimension (entities creating the products)
 *   the "source" dimension (entities being split down the middle)
 *   the "product" dimension (entities being created)
 * this function generates the "product" entities
 * which fill the open domain of the "element" entities
 * as a result of certain "source" entities being
 * split in the middle.
 *
 * to separate concerns, it relies on splits_to_elements
 * to preprocess some information about which elements
 * are splitting, which generated vertices to use when
 * making product entities, and the relative directions
 * from elements to source entities
 *
 * both elements and products are represented in terms
 * of their vertices
 */

void refine_topology(
    unsigned elem_dim,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    unsigned const* gen_offset_of_elems,
    unsigned const* gen_vert_of_elems,
    unsigned const* gen_direction_of_elems,
    unsigned* nprods_out,
    unsigned** verts_of_prods_out);

struct mesh;

void mesh_refine_topology(struct mesh* m,
    unsigned src_dim,
    unsigned prod_dim,
    unsigned const* gen_offset_of_elems,
    unsigned const* gen_vert_of_elems,
    unsigned const* gen_direction_of_elems,
    unsigned* nprods_out,
    unsigned** verts_of_prods_out);

#endif
