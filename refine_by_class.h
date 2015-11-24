#ifndef REFINE_BY_CLASS_H
#define REFINE_BY_CLASS_H

struct mesh;

/* Given a full topological representation
   with all entities classified, this
   function refines the mesh, putting new
   vertices in the center of old entities,
   until the output mesh is such that the
   classification of any mesh entity in
   it may be correctly inferred by the
   classification of its vertices,
   see infer_class().

   This is a somewhat hacky approach to meshing,
   and the alternative is to maintain a
   full representation.
*/

void refine_by_class(struct mesh** p_m);

#endif
