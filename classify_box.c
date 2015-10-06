#include "classify_box.h"
#include <math.h>    // for fabs
#include "loop.h"  // for malloc
#include "mesh.h"    // for mesh_add_nodal_label, mesh_count, mesh_find_noda...

#define EPSILON 1e-10

static void classify_box(
    unsigned nverts,
    double const* coords,
    unsigned** p_dims,
    unsigned** p_ids)
{
  unsigned* dims = loop_malloc(sizeof(unsigned) * nverts);
  unsigned* ids = loop_malloc(sizeof(unsigned) * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    double const* x = coords + i * 3;
    unsigned dim = 3;
    unsigned id = 0;
    for (unsigned fj = 0; fj < 3; ++fj) {
      unsigned j = 2 - fj;
      unsigned l = 1;
      if (fabs(x[j] - 0) < EPSILON)
        l = 0;
      else if (fabs(x[j] - 1) < EPSILON)
        l = 2;
      if (l != 1)
        --dim;
      id += l;
      id *= 3;
    }
    dims[i] = dim;
    ids[i] = id / 3;
  }
  *p_dims = dims;
  *p_ids = ids;
}

void mesh_classify_box(struct mesh* m)
{
  unsigned* dims;
  unsigned* ids;
  classify_box(mesh_count(m, 0),
    mesh_find_tag(m, 0, "coordinates")->data,
    &dims, &ids);
  mesh_add_tag(m, 0, TAG_U32, "class_dim", 1, dims);
  mesh_add_tag(m, 0, TAG_U32, "class_id", 1, ids);
}
