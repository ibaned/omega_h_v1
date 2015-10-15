#include "subset.h"

#include <string.h>

#include "ints.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "tables.h"
#include "tag.h"

#define GENERIC_SUBSET(T,name) \
T* name##_subset( \
    unsigned n, \
    unsigned width, \
    T const* a, \
    unsigned const* offsets) \
{ \
  unsigned nsub = offsets[n]; \
  T* out = LOOP_MALLOC(T, nsub * width); \
  for (unsigned i = 0; i < n; ++i) { \
    if (offsets[i] != offsets[i + 1]) \
      for (unsigned j = 0; j < width; ++j) \
        out[offsets[i] * width + j] = a[i * width + j]; \
  } \
  return out; \
}

GENERIC_SUBSET(unsigned char, uchars)
GENERIC_SUBSET(unsigned, uints)
GENERIC_SUBSET(unsigned long, ulongs)
GENERIC_SUBSET(double, doubles)

struct mesh* subset_mesh(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* offsets)
{
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned nelems_out = offsets[nelems];
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* verts_of_elems_out = uints_subset(nelems, verts_per_elem,
      verts_of_elems, offsets);
  unsigned nverts = mesh_count(m, 0);
  unsigned* marked_elems = uints_unscan(offsets, nelems);
  unsigned* marked_verts = mesh_mark_down(m, elem_dim, 0, marked_elems);
  loop_free(marked_elems);
  unsigned* vert_offsets = uints_exscan(marked_verts, nverts);
  loop_free(marked_verts);
  for (unsigned i = 0; i < nelems_out * verts_per_elem; ++i) {
    unsigned tmp = vert_offsets[verts_of_elems_out[i]];
    verts_of_elems_out[i] = tmp;
  }
  unsigned nverts_out = vert_offsets[nverts];
  struct mesh* out = new_mesh(elem_dim);
  mesh_set_ents(out, 0, nverts_out, 0);
  mesh_set_ents(out, elem_dim, nelems_out, verts_of_elems_out);
  for (unsigned i = 0; i < mesh_count_tags(m, 0); ++i) {
    struct const_tag* t = mesh_get_tag(m, 0, i);
    void* data_out = 0;
    switch (t->type) {
      case TAG_U8:  data_out = uchars_subset(
                        nverts, t->ncomps, t->d.u8, vert_offsets);
                    break;
      case TAG_U32: data_out = uints_subset(
                        nverts, t->ncomps, t->d.u32, vert_offsets);
                    break;
      case TAG_U64: data_out = ulongs_subset(
                        nverts, t->ncomps, t->d.u64, vert_offsets);
                    break;
      case TAG_F64: data_out = doubles_subset(
                        nverts, t->ncomps, t->d.f64, vert_offsets);
                    break;
    }
    mesh_add_tag(out, 0, t->type, t->name, t->ncomps, data_out);
  }
  loop_free(vert_offsets);
  return out;
}
