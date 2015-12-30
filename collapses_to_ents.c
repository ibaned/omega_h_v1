#include "collapses_to_ents.h"

#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

void collapses_to_ents(
    struct mesh* m,
    unsigned ent_dim,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned** p_gen_offset_of_ents,
    unsigned** p_gen_vert_of_ents,
    unsigned** p_gen_direction_of_ents,
    unsigned** p_offset_of_same_ents,
    unsigned** p_dead_sides)
{
  unsigned nents = mesh_count(m, ent_dim);
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned const* sides_of_ents = 0;
  unsigned* dead_sides = 0;
  if (p_dead_sides) {
    sides_of_ents = mesh_ask_down(m, ent_dim, ent_dim - 1);
    unsigned nsides = mesh_count(m, ent_dim - 1);
    dead_sides = uints_filled(nsides, 0);
  }
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned sides_per_ent = the_down_degrees[ent_dim][ent_dim - 1];
  unsigned* ent_will_gen = uints_filled(nents, 0);
  unsigned* gen_vert_of_ents = LOOP_MALLOC(unsigned, nents);
  unsigned* gen_direction_of_ents = LOOP_MALLOC(unsigned, nents);
  unsigned* ent_is_same = LOOP_MALLOC(unsigned, nents);
  unsigned const* ent_sides_opp_verts = the_opposite_orders[ent_dim][0];
  for (unsigned i = 0; i < nents; ++i) {
    unsigned const* verts_of_ent = verts_of_ents + i * verts_per_ent;
    unsigned col_vert = INVALID;
    unsigned direction = INVALID;
    for (unsigned j = 0; j < verts_per_ent; ++j) {
      unsigned vert = verts_of_ent[j];
      if (gen_offset_of_verts[vert] != gen_offset_of_verts[vert + 1]) {
        col_vert = vert;
        direction = j;
      }
    }
    if (direction == INVALID) {
      ent_is_same[i] = 1;
      /* not adjacent to collapsing vertex */
      continue;
    }
    ent_is_same[i] = 0;
    unsigned gen_vert = gen_vert_of_verts[col_vert];
    unsigned gen_dir = INVALID;
    for (unsigned j = 0; j < verts_per_ent; ++j)
      if (verts_of_ent[j] == gen_vert)
        gen_dir = j;
    unsigned ent_will_die = (gen_dir != INVALID);
    if (ent_will_die) {
      if (p_dead_sides) {
        unsigned ent_side = ent_sides_opp_verts[gen_dir];
        unsigned side = sides_of_ents[i * sides_per_ent + ent_side];
        dead_sides[side] = 1;
      }
      continue;
    }
    ent_will_gen[i] = 1;
    gen_vert_of_ents[i] = gen_vert;
    gen_direction_of_ents[i] = direction;
  }
  *p_gen_offset_of_ents = uints_exscan(ent_will_gen, nents);
  loop_free(ent_will_gen);
  *p_gen_vert_of_ents = gen_vert_of_ents;
  *p_gen_direction_of_ents = gen_direction_of_ents;
  *p_offset_of_same_ents = uints_exscan(ent_is_same, nents);
  loop_free(ent_is_same);
  if (p_dead_sides)
    *p_dead_sides = dead_sides;
}
