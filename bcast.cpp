#include "bcast.hpp"

#include <assert.h>
#include <string.h>

#include "comm.hpp"
#include "int_casts.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "tag.hpp"

static void bcast_tag(struct const_tag* t, struct tags* into)
{
  enum tag_type type = TAG_U8;
  if (!comm_rank())
    type = t->type;
  type = static_cast<enum tag_type>(comm_bcast_uint(U(type)));
  unsigned nl = 0;
  if (!comm_rank())
    nl = U(strlen(t->name));
  nl = comm_bcast_uint(nl);
  char* name = LOOP_HOST_MALLOC(char, nl + 1);
  if (!comm_rank())
    strcpy(name, t->name);
  comm_bcast_chars(name, nl + 1);
  unsigned ncomps = 0;
  if (!comm_rank())
    ncomps = t->ncomps;
  ncomps = comm_bcast_uint(ncomps);
  if (comm_rank())
    add_tag(into, type, name, ncomps, LOOP_MALLOC(char, 0));
  loop_host_free(name);
}

static void bcast_tags(struct tags* from, struct tags* into)
{
  unsigned ntags = 0;
  if (!comm_rank())
    ntags = count_tags(from);
  ntags = comm_bcast_uint(ntags);
  for (unsigned i = 0; i < ntags; ++i) {
    struct const_tag* t = 0;
    if (!comm_rank())
      t = get_tag(from, i);
    bcast_tag(t, into);
  }
}

struct mesh* bcast_mesh_metadata(struct mesh* m)
{
  if (!comm_rank())
    assert(m);
  else
    assert(!m);
  unsigned dim = 0;
  if (!comm_rank())
    dim = mesh_dim(m);
  dim = comm_bcast_uint(dim);
  enum mesh_rep rep = MESH_REDUCED;
  if (!comm_rank())
    rep = mesh_get_rep(m);
  rep = static_cast<enum mesh_rep>(comm_bcast_uint(U(rep)));
  if (comm_rank())
    m = new_mesh(dim, rep, 0);
  for (unsigned i = 0; i <= dim; ++i) {
    unsigned has = 0;
    if (!comm_rank())
      has = mesh_has_dim(m, i);
    has = comm_bcast_uint(has);
    if (!has)
      continue;
    if (comm_rank())
      mesh_set_ents(m, i, 0, LOOP_MALLOC(unsigned, 0));
    struct tags* ts = mesh_tags(m, i);
    bcast_tags(ts, ts);
  }
  return m;
}
