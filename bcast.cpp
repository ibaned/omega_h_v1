#include "bcast.hpp"

#include <cassert>
#include <cstring>

#include "comm.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tag.hpp"
#include "int_casts.hpp"

namespace omega_h {

static void bcast_tag(struct const_tag* t, struct tags* into, unsigned is_source)
{
  enum tag_type type = TAG_U8;
  if (is_source)
    type = t->type;
  type = static_cast<enum tag_type>(comm_bcast_uint(U(type)));
  unsigned nl = 0;
  if (is_source)
    nl = U(strlen(t->name));
  nl = comm_bcast_uint(nl);
  char* name = LOOP_HOST_MALLOC(char, nl + 1);
  if (is_source)
    strcpy(name, t->name);
  comm_bcast_chars(name, nl + 1);
  unsigned ncomps = 0;
  if (is_source)
    ncomps = t->ncomps;
  ncomps = comm_bcast_uint(ncomps);
  if (!is_source)
    add_tag(into, type, name, ncomps, LOOP_MALLOC(char, 0));
  loop_host_free(name);
}

static void bcast_tags(struct tags* from, struct tags* into, unsigned is_source)
{
  unsigned ntags = 0;
  if (is_source)
    ntags = count_tags(from);
  ntags = comm_bcast_uint(ntags);
  for (unsigned i = 0; i < ntags; ++i) {
    struct const_tag* t = 0;
    if (is_source)
      t = get_tag(from, i);
    bcast_tag(t, into, is_source);
  }
}

struct mesh* bcast_mesh_metadata(struct mesh* m, unsigned is_source)
{
  if (is_source)
    assert(m);
  else
    assert(!m);
  unsigned dim = 0;
  if (is_source)
    dim = mesh_dim(m);
  dim = comm_bcast_uint(dim);
  enum mesh_rep rep = MESH_REDUCED;
  if (is_source)
    rep = mesh_get_rep(m);
  rep = static_cast<enum mesh_rep>(comm_bcast_uint(U(rep)));
  if (!is_source)
    m = new_mesh(dim, rep, 0);
  for (unsigned i = 0; i <= dim; ++i) {
    unsigned has = 0;
    if (is_source)
      has = mesh_has_dim(m, i);
    has = comm_bcast_uint(has);
    if (!has)
      continue;
    if (!is_source)
      mesh_set_ents(m, i, 0, LOOP_MALLOC(unsigned, 0));
    struct tags* ts = mesh_tags(m, i);
    bcast_tags(ts, ts, is_source);
  }
  return m;
}

}
