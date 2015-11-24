#include "refine_by_class.h"

#include <assert.h>

#include "copy_mesh.h"
#include "find_by_verts.h"
#include "infer_class.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "refine_common.h"

/* Multiple assumptions come together here.
   First, we are using the dimension-generic
   entity splitting system in refine_common().
   each candidate entity will be split in the
   barycenter by a new vertex.

   we will only mark candidates whose inferred
   and real classifications do no not match.
   once refined, the open domain of that
   entity will be resolved, since all the refined
   entities that fill it will be correctly
   classified by the new vertex.
   this relies on refine_class() to accept existing
   classification arrays.

   we will proceed from high to low dimension,
   first splitting elements whose classification
   differs from their vertices, and edges will
   be split last. the reason for this is that
   splitting an entity in the middle does not
   affect its boundary, and therefore we can be
   sure that all entities will be properly corrected.

   since refine_common() only removes entities
   greater than or equal to the source dimension,
   we can assume that the (n) dimensional simplices
   from the input mesh will exist
   until we refine (n) dimensional sources.
   they may have be re-derived however, so we search
   by corner vertices for matches.
*/

void refine_by_class(struct mesh** p_m)
{
  /* copy the original mesh */
  struct mesh* m_in = copy_mesh(*p_m);
  unsigned dim = mesh_dim(m_in);
  for (unsigned rd = 0; rd < dim; ++rd) {
  /* from elements to edges, not vertices */
    unsigned d = dim - rd;
    unsigned nfixed = 0;
    while (1) {
      unsigned nents_in = mesh_count(m_in, d);
      unsigned nents = mesh_count(*p_m, d);
      unsigned* candidates = uints_filled(nents, 0);
      /* force inferrence of current classifications and make copies for
         modification */
      if (mesh_find_tag(*p_m, d, "class_dim"))
        mesh_free_tag(*p_m, d, "class_dim");
      if (mesh_find_tag(*p_m, d, "class_id"))
        mesh_free_tag(*p_m, d, "class_id");
      unsigned* class_dims = uints_copy(mesh_ask_class_dim(*p_m, d), nents);
      unsigned* class_ids = uints_copy(mesh_ask_class_id(*p_m, d), nents);
      /* get the original entities and their classifications */
      unsigned const* verts_of_ents_in = mesh_ask_down(m_in, d, 0);
      unsigned const* class_dims_in = mesh_find_tag(m_in, d, "class_dim")->d.u32;
      unsigned const* class_ids_in = mesh_find_tag(m_in, d, "class_id")->d.u32;
      unsigned verts_per_ent = the_down_degrees[d][0];
      unsigned const* verts_of_ents = mesh_ask_down(*p_m, d, 0);
      unsigned const* ents_of_verts = mesh_ask_up(*p_m, 0, d)->adj;
      unsigned const* ents_of_verts_offsets = mesh_ask_up(*p_m, 0, d)->offsets;
      unsigned nskipped = 0;
      unsigned ncandidates = 0;
      for (unsigned ent_in = 0; ent_in < nents_in; ++ent_in) {
        unsigned ent = find_by_verts(verts_per_ent,
            verts_of_ents_in + ent_in * verts_per_ent,
            verts_of_ents, ents_of_verts, ents_of_verts_offsets);
        if (ent == INVALID) {
          ++nskipped;
          continue;
        }
        candidates[ent] = ((class_dims[ent] != class_dims_in[ent_in]) ||
                           (class_ids[ent] != class_ids_in[ent_in]));
        if (candidates[ent])
          ++ncandidates;
        class_dims[ent] = class_dims_in[ent_in];
        class_ids[ent] = class_ids_in[ent_in];
      }
      assert(nskipped == nfixed);
      /* free the current classifications and replace them
         with the modified classifications, which will
         induce refinement to correctly classify new vertices */
      mesh_free_tag(*p_m, d, "class_dim");
      mesh_add_tag(*p_m, d, TAG_U32, "class_dim", 1, class_dims);
      mesh_free_tag(*p_m, d, "class_id");
      mesh_add_tag(*p_m, d, TAG_U32, "class_id", 1, class_ids);
      unsigned nverts_before = mesh_count(*p_m, 0);
      /* refine without restrictions. we need inferrable classification,
         possibly at heavy cost to element quality... */
      refine_common(p_m, d, candidates, 0.0, 0);
      loop_free(candidates);
      unsigned nverts_after = mesh_count(*p_m, 0);
      nfixed += nverts_after - nverts_before;
      if (nverts_after - nverts_before == ncandidates)
        break;
    }
  }
  free_mesh(m_in);
}
