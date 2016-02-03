#include "collapses_to_ents.h"

#include "arrays.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

LOOP_KERNEL(collapse_to_ent,
    unsigned const* verts_of_ents,
    unsigned verts_per_ent,
    unsigned sides_per_ent,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned* ent_sides_opp_verts,
    unsigned const* sides_of_ents,
    unsigned const* fused_ents,
    unsigned* ent_is_same,
    unsigned* fused_sides,
    unsigned* ent_will_gen,
    unsigned* gen_vert_of_ents,
    unsigned* gen_direction_of_ents)

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
    return;
  }
  ent_is_same[i] = 0;
  unsigned gen_vert = gen_vert_of_verts[col_vert];
  unsigned gen_dir = INVALID;
  for (unsigned j = 0; j < verts_per_ent; ++j)
    if (verts_of_ent[j] == gen_vert)
      gen_dir = j;
  unsigned ent_will_die = (gen_dir != INVALID);
  if (ent_will_die) {
    if (fused_sides) {
      unsigned ent_side = ent_sides_opp_verts[gen_dir];
      unsigned side = sides_of_ents[i * sides_per_ent + ent_side];
      fused_sides[side] = 1;
    }
    return;
  }
  unsigned is_fused_side = (fused_ents && fused_ents[i]);
  if (!is_fused_side) {
    ent_will_gen[i] = 1;
    gen_vert_of_ents[i] = gen_vert;
    gen_direction_of_ents[i] = direction;
  }
}

void collapses_to_ents(
    struct mesh* m,
    unsigned ent_dim,
    unsigned const* gen_offset_of_verts,
    unsigned const* gen_vert_of_verts,
    unsigned const* fused_ents,
    unsigned** p_gen_offset_of_ents,
    unsigned** p_gen_vert_of_ents,
    unsigned** p_gen_direction_of_ents,
    unsigned** p_offset_of_same_ents,
    unsigned** p_fused_sides)
{
  unsigned nents = mesh_count(m, ent_dim);
  unsigned const* verts_of_ents = mesh_ask_down(m, ent_dim, 0);
  unsigned const* sides_of_ents = 0;
  unsigned* fused_sides = 0;
  if (p_fused_sides) {
    sides_of_ents = mesh_ask_down(m, ent_dim, ent_dim - 1);
    unsigned nsides = mesh_count(m, ent_dim - 1);
    fused_sides = uints_filled(nsides, 0);
  }
  unsigned verts_per_ent = the_down_degrees[ent_dim][0];
  unsigned sides_per_ent = the_down_degrees[ent_dim][ent_dim - 1];
  unsigned* ent_will_gen = uints_filled(nents, 0);
  unsigned* gen_vert_of_ents = LOOP_MALLOC(unsigned, nents);
  unsigned* gen_direction_of_ents = LOOP_MALLOC(unsigned, nents);
  unsigned* ent_is_same = LOOP_MALLOC(unsigned, nents);
  unsigned* ent_sides_opp_verts = LOOP_TO_DEVICE(unsigned,
      the_opposite_orders[ent_dim][0], verts_per_ent);
  LOOP_EXEC(collapse_to_ent, nents,
      verts_of_ents,
      verts_per_ent,
      sides_per_ent,
      gen_offset_of_verts,
      gen_vert_of_verts,
      ent_sides_opp_verts,
      sides_of_ents,
      fused_ents,
      ent_is_same,
      fused_sides,
      ent_will_gen,
      gen_vert_of_ents,
      gen_direction_of_ents);
  loop_free(ent_sides_opp_verts);
  *p_gen_offset_of_ents = uints_exscan(ent_will_gen, nents);
  loop_free(ent_will_gen);
  *p_gen_vert_of_ents = gen_vert_of_ents;
  *p_gen_direction_of_ents = gen_direction_of_ents;
  *p_offset_of_same_ents = uints_exscan(ent_is_same, nents);
  loop_free(ent_is_same);
  if (p_fused_sides)
    *p_fused_sides = fused_sides;
}
