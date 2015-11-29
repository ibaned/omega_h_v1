#include "migrate_cloud.h"

#include <assert.h>

#include "cloud.h"
#include "exchanger.h"
#include "loop.h"

#define GENERIC_MIGRATE(T, name) \
static T* name##_migrate(struct exchanger* ex, unsigned ncomps, \
    T const* data) \
{ \
  T* recvd = LOOP_MALLOC(T, ex->nrecvd * ncomps); \
  for (unsigned i = 0; i < ex->ndests; ++i) { \
    assert(ex->recvd_of_dests_offsets[i] == i); \
    unsigned irecvd = ex->recvd_of_dests[i]; \
    for (unsigned j = 0; j < ncomps; ++j) \
      recvd[irecvd * ncomps + j] = data[i * ncomps + j]; \
  } \
  T* out = unexchange_##name(ex, ncomps, recvd); \
  loop_free(recvd); \
  return out; \
}

GENERIC_MIGRATE(unsigned, uints)
GENERIC_MIGRATE(unsigned long, ulongs)
GENERIC_MIGRATE(double, doubles)

void migrate_cloud(struct cloud** p_c,
    unsigned nrecvd,
    unsigned const* recvd_ranks,
    unsigned const* recvd_ids)
{
  struct cloud* c = *p_c;
  struct exchanger* ex = new_exchanger(
      nrecvd, cloud_count(c),
      recvd_ranks, recvd_ids);
  struct cloud* c_out = new_cloud(nrecvd);
  for (unsigned i = 0; i < cloud_count_tags(c); ++i) {
    struct const_tag* t = cloud_get_tag(c, i);
    void* data_out = 0;
    switch (t->type) {
      case TAG_U8:
        data_out = 0;
        break;
      case TAG_U32:
        data_out = uints_migrate(ex, t->ncomps, t->d.u32);
        break;
      case TAG_U64:
        data_out = ulongs_migrate(ex, t->ncomps, t->d.u64);
        break;
      case TAG_F64:
        data_out = doubles_migrate(ex, t->ncomps, t->d.f64);
        break;
    }
    assert(data_out);
    cloud_add_tag(c_out, t->type, t->name, t->ncomps, data_out);
  }
  free_exchanger(ex);
  free_cloud(*p_c);
  *p_c = c_out;
}
