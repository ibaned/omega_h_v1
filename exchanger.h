#ifndef EXCHANGER_H
#define EXCHANGER_H

struct comm;

/* this system encompasses the exchange of many small items
   to many small destinations on other MPI ranks.
   for example, consider sending values from duplicate nodes
   of a mesh to the "owner" or "master" nodes
*/

struct exchanger {
  struct comm* forward_comm;
  struct comm* reverse_comm;
  unsigned nsent;
  unsigned nrecvd;
  unsigned nsends;
  unsigned nrecvs;
  unsigned ndests;
  unsigned padding__;
  unsigned* send_ranks;
  unsigned* recv_ranks;
  unsigned* send_counts;
  unsigned* recv_counts;
  unsigned* send_offsets;
  unsigned* recv_offsets;
  unsigned* send_of_sent;
  unsigned* send_idx_of_sent;
  unsigned* recv_of_recvd;
  unsigned* recvd_of_dests;
  unsigned* recvd_of_dests_offsets;
};

struct exchanger* new_exchanger(
    unsigned nsent, /* number of items to send */
    /* number of destinations on this MPI rank
       for items from other MPI ranks */
    unsigned ndests,
    /* destinations of each item */
    unsigned const* dest_rank_of_sent,
    unsigned const* dest_idx_of_sent);
unsigned* exchange_uints(struct exchanger* ex, unsigned width,
    unsigned const* sent);
double* exchange_doubles(struct exchanger* ex, unsigned width,
    double const* sent);
unsigned* unexchange_uints(struct exchanger* ex, unsigned width,
    unsigned const* recvd);
double* unexchange_doubles(struct exchanger* ex, unsigned width,
    double const* recvd);
unsigned long* unexchange_ulongs(struct exchanger* ex, unsigned width,
    unsigned long const* recvd);
void free_exchanger(struct exchanger* ex);

#endif
