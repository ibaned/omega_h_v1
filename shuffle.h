#ifndef SHUFFLE_H
#define SHUFFLE_H

struct shuffle;

struct shuffle* new_shuffle(unsigned n,
    unsigned const* parts, unsigned const* indices);
void print_shuffle(struct shuffle* s);
void free_shuffle(struct shuffle* s);

unsigned shuffle_out_size(struct shuffle* s);
void shuffle_uints(struct shuffle* s, unsigned const* in, unsigned** out);
unsigned const* shuffle_offsets(struct shuffle* s);

#endif
