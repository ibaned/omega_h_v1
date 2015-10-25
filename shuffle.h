#ifndef SHUFFLE_H
#define SHUFFLE_H

struct shuffle;

struct shuffle* new_shuffle(unsigned n, unsigned const* parts);
void print_shuffle(struct shuffle* s);
void free_shuffle(struct shuffle* s);
unsigned* shuffle_uints(struct shuffle* s, unsigned const* a);
unsigned shuffle_nrecv_ents(struct shuffle* s);

#endif
