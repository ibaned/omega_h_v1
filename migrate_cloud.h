#ifndef MIGRATE_CLOUD_H
#define MIGRATE_CLOUD_H

/* pull-based migration of a cloud.
   each rank indicates which particles
   from the input cloud it would like
   to have in the output cloud */

struct cloud;

void migrate_cloud(struct cloud** p_c,
    unsigned nrecvd,
    unsigned const* recvd_ranks,
    unsigned const* recvd_ids);

#endif
