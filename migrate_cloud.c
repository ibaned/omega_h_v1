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
  struct exchanger* push = make_reverse_exchanger(cloud_count(c), nrecvd,
      recvd_ranks, recvd_ids);
  struct cloud* c_out = new_cloud(nrecvd);
  push_tags(push, cloud_tags(c), cloud_tags(c_out));
  free_exchanger(push);
  free_cloud(*p_c);
  *p_c = c_out;
}
