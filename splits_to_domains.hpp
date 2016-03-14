#ifndef SPLITS_TO_ELEMENTS_HPP
#define SPLITS_TO_ELEMENTS_HPP

namespace omega_h {

/* this is a preprocessing function for
 * refine_topology.
 * the entities whose interior domains are being
 * filled are referred to as domain entities (dom_dim),
 * and the entities that are getting new vertices placed
 * in their centers (splitting) are source entities.
 * each domain entity examines its adjacent source
 * entities.
 * at most one of these will actually be splitting,
 * indicated by offset_of_srcs.
 * domain entities adjacent to splitting sources
 * are marked (offset_of_doms), and the relative
 * direction and new vertex of the adjacent source
 * are recorded (direction_of_doms, vert_of_doms).
 */

void project_splits_to_domains(
    unsigned dom_dim,
    unsigned src_dim,
    unsigned nelems,
    unsigned const* srcs_of_doms,
    unsigned const* offset_of_srcs,
    unsigned const* vert_of_srcs,
    unsigned** p_offset_of_doms,
    unsigned** p_direction_of_doms,
    unsigned** p_vert_of_doms);

struct mesh;

void mesh_splits_to_domains(
  struct mesh* m,
  unsigned dom_dim,
  unsigned src_dim,
  unsigned const* offset_of_srcs,
  unsigned const* vert_of_srcs,
  unsigned** p_offset_of_doms,
  unsigned** p_direction_of_doms,
  unsigned** p_vert_of_doms);

}

#endif
