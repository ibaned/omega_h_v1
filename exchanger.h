#ifndef EXCHANGER_H
#define EXCHANGER_H

struct comm;

struct exchanger {
  struct comm* forward_comm;
  struct comm* reverse_comm;
  unsigned nsent;
  unsigned nrecvd;
  unsigned nsends;
  unsigned nrecvs;
  unsigned* send_parts;
  unsigned* recv_parts;
  unsigned* send_counts;
  unsigned* recv_counts;
  unsigned* send_offsets;
  unsigned* recv_offsets;
  unsigned* send_of_sent;
  unsigned* send_idx_of_sent;
};

struct exchanger* new_exchanger(unsigned nsent, unsigned const* part_of_sent);
unsigned* exchange_uints(struct exchanger* ex, unsigned width,
    unsigned const* sent);
unsigned* unexchange_uints(struct exchanger* ex, unsigned width,
    unsigned const* recvd);
void free_exchanger(struct exchanger* ex);

#endif
