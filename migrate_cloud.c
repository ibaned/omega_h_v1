#include "migrate_cloud.h"

#include <assert.h>

#include "cloud.h"
#include "exchanger.h"
#include "loop.h"

void migrate_cloud(struct cloud** p_c,
    unsigned nrecvd,
    unsigned const* recvd_ranks,
    unsigned const* recvd_ids)
{
  struct cloud* c = *p_c;
  struct exchanger* ex = new_exchanger(nrecvd, recvd_ranks);
  set_exchanger_dests(ex, cloud_count(c), recvd_ids);
  struct cloud* c_out = new_cloud(nrecvd);
  for (unsigned i = 0; i < cloud_count_tags(c); ++i) {
    struct const_tag* t = cloud_get_tag(c, i);
    void* data_out = 0;
    switch (t->type) {
      case TAG_U8:
        data_out = 0;
        break;
      case TAG_U32:
        data_out = exchange_uints(ex, t->ncomps, t->d.u32,
            EX_REV, EX_ITEM);
        break;
      case TAG_U64:
        data_out = exchange_ulongs(ex, t->ncomps, t->d.u64,
            EX_REV, EX_ITEM);
        break;
      case TAG_F64:
        data_out = exchange_doubles(ex, t->ncomps, t->d.f64,
            EX_REV, EX_ITEM);
        break;
    }
    assert(data_out);
    cloud_add_tag(c_out, t->type, t->name, t->ncomps, data_out);
  }
  free_exchanger(ex);
  free_cloud(*p_c);
  *p_c = c_out;
}
