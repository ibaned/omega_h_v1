#include "eval_field.h"
#include "loop.h"  // for malloc
#include "field.h"   // for const_field
#include "mesh.h"    // for mesh_add_nodal_field, mesh_count, mesh_find_noda...

double* eval_field(
    unsigned nverts,
    double const* coords,
    unsigned ncomps,
    void (*fun)(double const x[3], double out[]))
{
  double* out = loop_malloc(sizeof(double) * ncomps * nverts);
  double const* vert_coords = coords;
  double* vert_comps = out;
  for (unsigned i = 0; i < nverts; ++i) {
    fun(vert_coords, vert_comps);
    vert_coords += 3;
    vert_comps += ncomps;
  }
  return out;
}

void mesh_eval_field(struct mesh* m, char const* name, unsigned ncomps,
    void (*fun)(double const x[3], double out[]))
{
  double* data = eval_field(mesh_count(m, 0),
      mesh_find_field(m, 0, "coordinates")->data, ncomps, fun);
  mesh_add_field(m, 0, name, ncomps, data);
}
