#include "eval_field.hpp"

#include "arrays.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tag.hpp"

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
  double* host_coords = doubles_to_host(coords, nents * 3);
  double* host_out = LOOP_HOST_MALLOC(double, ncomps * nents);
  for (unsigned i = 0; i < nents; ++i) {
    double const* ent_coords = host_coords + i * 3;
    double* ent_out = host_out + i * ncomps;
    fun(ent_coords, ent_out);
  }
  loop_host_free(host_coords);
  double* out = array_to_device(host_out, ncomps * nents);
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

void mesh_eval_field2(struct mesh* m, unsigned ent_dim, char const* name,
    unsigned ncomps, enum osh_transfer tt,
    void (*fun)(double const* x, double* out))
{
  double* data = eval_field(mesh_count(m, ent_dim),
      mesh_find_tag(m, ent_dim, "coordinates")->d.f64, ncomps, fun);
  add_tag2(mesh_tags(m, ent_dim), TAG_F64, name, ncomps, tt, data);
}
