#ifndef EXCHANGER_H
#define EXCHANGER_H

struct comm;

enum exch_dir {
  EX_FOR,
  EX_REV,
};

enum exch_start {
  EX_ROOT,
  EX_ITEM,
};

#define EX_DIRS 2

/* The exchanger is designed to send information
   from many small sources to many small destinations.
   It is a symmetric multi-stage pipeline as follows.

   1) expand source array to input array.
      the offsets should specify a one-to-many
      mapping from source items to input items.

   2) shuffle the input array such that items are
      sorted by destination rank
      (the result is the sent array)

   3) use comm_exch_* to exchange the sent array data,
      resulting in a source-rank-sorted received array

   4) shuffle the received array, resulting in the
      "output" array which is sorted by destination
      item

   5) maintain the one-to-many mapping between destination
      items and output items, this is just offered to
      the user for higher-level purposes

   The above five stages can be reversed, thus the symmetric
   nature of the exchanger.
   For now, the exchanger is constructed in a push-based
   fashion, by giving, for each *input* item, the
   global ID (rank and local id) of the *destination* item,
   and optionally the offsets defining the
   source to input mapping.
*/

struct exchanger {
  struct comm* comms[EX_DIRS]; /* forward and reverse communicators */
  unsigned nitems[EX_DIRS]; /* sent and received item counts */
  unsigned nroots[EX_DIRS]; /* source and destination counts */
  unsigned nmsgs[EX_DIRS]; /* number of unique neighbor ranks */
  unsigned* ranks[EX_DIRS]; /* neighbor ranks */
  unsigned* msg_counts[EX_DIRS]; /* message sizes */
/* offsets into message-sorted data arrays of each message */
  unsigned* msg_offsets[EX_DIRS];
/* 1-to-1 maps from root-sorted to message-sorted data */
  unsigned* shuffles[EX_DIRS];
/* message index for each input and output item */
  unsigned* msg_of_items[EX_DIRS];
/* one-to-many mapping from roots to items */
  unsigned* items_of_roots_offsets[EX_DIRS];
};

struct exchanger* new_exchanger(
    unsigned nsent,
    unsigned const* dest_rank_of_sent);

void set_exchanger_dests(
    struct exchanger* ex,
    /* number of destinations on this MPI rank
       for items from other MPI ranks */
    unsigned ndests,
    unsigned const* dest_idx_of_sent);

void set_exchanger_srcs(
    struct exchanger* ex,
    unsigned nsrcs,
    /* map from sources to input items */
    unsigned const* sent_of_srcs_offsets);

unsigned* exchange_uints(struct exchanger* ex, unsigned width,
    unsigned const* data, enum exch_dir dir, enum exch_start start);
double* exchange_doubles(struct exchanger* ex, unsigned width,
    double const* data, enum exch_dir dir, enum exch_start start);
unsigned long* exchange_ulongs(struct exchanger* ex, unsigned width,
    unsigned long const* data, enum exch_dir dir, enum exch_start start);

void free_exchanger(struct exchanger* ex);

void reverse_exchanger(struct exchanger* ex);

struct exchanger* make_reverse_exchanger(unsigned nsent, unsigned nrecvd,
    unsigned const* recvd_ranks, unsigned const* recvd_ids);

#endif
