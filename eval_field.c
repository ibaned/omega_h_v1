#include "eval_field.h"
#include "loop.h"  // for malloc
#include "mesh.h"    // for mesh_add_nodal_field, mesh_count, mesh_find_noda...

double* eval_field(
    unsigned nents,
    double const* coords,
    unsigned ncomps,
    void (*fun)(double const x[3], double out[]))
{
  double* out = loop_malloc(sizeof(double) * ncomps * nents);
  for (unsigned i = 0; i < nents; ++i) {
    double const* ent_coords = coords + i * 3;
    double* ent_out = out + i * ncomps;
    fun(ent_coords, ent_out);
  }
  return out;
}

void mesh_eval_field(struct mesh* m, unsigned ent_dim, char const* name,
    unsigned ncomps, void (*fun)(double const x[3], double out[]))
{
  double* data = eval_field(mesh_count(m, ent_dim),
      mesh_find_tag(m, ent_dim, "coordinates")->data, ncomps, fun);
  mesh_add_tag(m, ent_dim, TAG_F64, name, ncomps, data);
}
