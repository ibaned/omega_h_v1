#ifndef COLLAPSES_TO_ENTS_H
#define COLLAPSES_TO_ENTS_H

struct mesh;

/* converts a vertex-wise description of
   what collapses to perform into an entity-wise
   description. */

void collapses_to_ents(
    struct mesh* m,
/* dimension of entities to project onto */
    unsigned ent_dim,
/* offsets denoting vertices being removed by
   edge collapses */
    unsigned const* gen_offset_of_verts,
/* for vertices removed by an edge collapse,
   maps them to the vertex they collapse onto */
    unsigned const* gen_vert_of_verts,
/* indicates whether an entity of (ent_dim) is "fused"
   in a collapse, there are several entities that are
   "collapsing onto" other entities, removing the entity
   between them.
   the collapsing vertex is the first; it collapses onto
   another vertex and removes the edge in between.
   now consider a triangle adjacent to the collapsing edge.
   one of the other two edges of that triangle is being
   collapsed onto the other edge, "pinching" the triangle
   to zero area.
   that edge (the one adjacent to the collapsing vertex)
   is a "fused" entity.
   these entities are a special category because they are
   adjacent to the collapsing vertex but they
   should *not* generate new entities because the resulting
   entity already exists as part of the boundary. */
    unsigned const* fused_ents,
/* offsets describing old entities that will generate new entities */
    unsigned** gen_offset_of_ents_out,
/* for each generating entity, maps it to the new vertex */
    unsigned** gen_vert_of_ents_out,
/* for generating entities, which vertex is it that is moving */
    unsigned** gen_direction_of_ents_out,
/* offsets describing entities that will remain in the output mesh */
    unsigned** offset_of_same_ents_out,
/* this function also does the job of identifying the "fused" entities
   of dimension (ent_dim - 1), called fused sides.
   this is later fed back to this function as the "fused_ents" argument,
   requiring that the function be called in decreasing order of (ent_dim) */
    unsigned** p_fused_sides);

#endif
