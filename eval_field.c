#include "eval_field.h"

#include "arrays.h"
#include "cloud.h"
#include "loop.h"
#include "mesh.h"
#include "tag.h"

/* for now, we carry this out on the host.
   if we want this done on the device, then
   we'll have to use a __device__ function
   pointer and probably relocatable device code */

double* eval_field(
    unsigned nents,
    double const* coords,
    unsigned ncomps,
    void (*fun)(double const* x, double* out))
{
  double* host_coords = LOOP_TO_HOST(double, coords, nents * 3);
  double* host_out = LOOP_HOST_MALLOC(double, ncomps * nents);
  for (unsigned i = 0; i < nents; ++i) {
    double const* ent_coords = host_coords + i * 3;
    double* ent_out = host_out + i * ncomps;
    fun(ent_coords, ent_out);
  }
  loop_host_free(host_coords);
  double* out = doubles_to_device(host_out, ncomps * nents);
  loop_host_free(host_out);
  return out;
}

void mesh_eval_field(struct mesh* m, unsigned ent_dim, char const* name,
    unsigned ncomps, void (*fun)(double const* x, double* out))
{
  double* data = eval_field(mesh_count(m, ent_dim),
      mesh_find_tag(m, ent_dim, "coordinates")->d.f64, ncomps, fun);
  mesh_add_tag(m, ent_dim, TAG_F64, name, ncomps, data);
}

void cloud_eval_field(struct cloud* m, char const* name,
    unsigned ncomps, void (*fun)(double const* x, double* out))
{
  double* data = eval_field(cloud_count(m),
      cloud_find_tag(m, "coordinates")->d.f64, ncomps, fun);
  cloud_add_tag(m, TAG_F64, name, ncomps, data);
}
