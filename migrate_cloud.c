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
  exchange_tags(ex, cloud_tags(c), cloud_tags(c_out), EX_REV, EX_ITEM);
  free_exchanger(ex);
  free_cloud(*p_c);
  *p_c = c_out;
}
