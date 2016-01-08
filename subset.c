#include "subset.h"

#include <string.h>

#include "arrays.h"
#include "ints.h"
#include "loop.h"
#include "mark.h"
#include "mesh.h"
#include "tables.h"
#include "tag.h"

void tags_subset(struct mesh* in, struct mesh* out,
    unsigned dim, unsigned const* offsets)
{
  unsigned nverts = mesh_count(in, dim);
  for (unsigned i = 0; i < mesh_count_tags(in, dim); ++i) {
    struct const_tag* t = mesh_get_tag(in, dim, i);
    void* vals_out = 0;
    switch (t->type) {
      case TAG_U8:
        vals_out = uchars_expand(nverts, t->ncomps, t->d.u8,
            offsets);
        break;
      case TAG_U32:
        vals_out = uints_expand(nverts, t->ncomps, t->d.u32,
            offsets);
        break;
      case TAG_U64:
        vals_out = ulongs_expand(nverts, t->ncomps, t->d.u64,
            offsets);
        break;
      case TAG_F64:
        vals_out = doubles_expand(nverts, t->ncomps, t->d.f64,
            offsets);
        break;
    }
    mesh_add_tag(out, dim, t->type, t->name, t->ncomps, vals_out);
  }
}

struct mesh* subset_mesh(
    struct mesh* m,
    unsigned elem_dim,
    unsigned const* offsets)
{
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned nelems_out = offsets[nelems];
  unsigned const* verts_of_elems = mesh_ask_down(m, elem_dim, 0);
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* verts_of_elems_out = uints_expand(nelems, verts_per_elem,
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
  /* TODO: should subset_mesh preserve full topology ? */
  struct mesh* out = new_mesh(elem_dim, MESH_REDUCED);
  mesh_set_ents(out, 0, nverts_out, 0);
  mesh_set_ents(out, elem_dim, nelems_out, verts_of_elems_out);
  tags_subset(m, out, 0, vert_offsets);
  tags_subset(m, out, elem_dim, offsets);
  loop_free(vert_offsets);
  return out;
}

void subset_verts_of_doms(
    struct mesh* m,
    unsigned dom_dim,
    unsigned const* offset_of_doms,
    unsigned* verts_of_prods)
{
  unsigned const* verts_of_doms = mesh_ask_down(m, dom_dim, 0);
  unsigned* offset_of_same = uints_negate_offsets(offset_of_doms,
      mesh_count(m, dom_dim));
  uints_expand_into(mesh_count(m, dom_dim), the_down_degrees[dom_dim][0],
      verts_of_doms, offset_of_same, verts_of_prods);
  loop_free(offset_of_same);
}
