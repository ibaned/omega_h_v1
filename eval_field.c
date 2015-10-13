#include "eval_field.h"

#include "cloud.h"
#include "loop.h"
#include "mesh.h"
#include "tag.h"

double* eval_field(
    unsigned nents,
    double const* coords,
    unsigned ncomps,
    void (*fun)(double const* x, double* out))
{
  double* out = LOOP_MALLOC(double, ncomps * nents);
  for (unsigned i = 0; i < nents; ++i) {
    double const* ent_coords = coords + i * 3;
    double* ent_out = out + i * ncomps;
    fun(ent_coords, ent_out);
  }
  return out;
}

void mesh_eval_field(struct mesh* m, unsigned ent_dim, char const* name,
    unsigned ncomps, void (*fun)(double const* x, double* out))
{
  double* data = eval_field(mesh_count(m, ent_dim),
      mesh_find_tag(m, ent_dim, "coordinates")->data, ncomps, fun);
  mesh_add_tag(m, ent_dim, TAG_F64, name, ncomps, data);
}

void cloud_eval_field(struct cloud* m, char const* name,
    unsigned ncomps, void (*fun)(double const* x, double* out))
{
  double* data = eval_field(cloud_count(m),
      cloud_find_tag(m, "coordinates")->data, ncomps, fun);
  cloud_add_tag(m, TAG_F64, name, ncomps, data);
}
