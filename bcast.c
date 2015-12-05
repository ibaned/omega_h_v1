#include "bcast.h"

#include <assert.h>
#include <string.h>

#include "comm.h"
#include "loop.h"
#include "mesh.h"
#include "tag.h"

static void bcast_tag(struct const_tag* t, struct tags* into)
{
  enum tag_type type = 0;
  if (!comm_rank())
    type = t->type;
  type = (enum tag_type) (comm_bcast_uint((unsigned) type));
  unsigned nl = 0;
  if (!comm_rank())
    nl = (unsigned) strlen(t->name);
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
  for (unsigned i = 0; i < ntags; ++i)
    bcast_tag(get_tag(from, i), into);
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
  if (comm_rank())
    m = new_mesh(dim);
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
