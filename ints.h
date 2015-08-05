#ifndef INTS_H
#define INTS_H

void ints_zero(unsigned* a, unsigned n);
void ints_copy(unsigned const* a, unsigned* b, unsigned n);
unsigned ints_max(unsigned const* a, unsigned n);
unsigned* ints_exscan(unsigned const* a, unsigned n);
unsigned* ints_unscan(unsigned const* a, unsigned n);
unsigned* ints_negate(unsigned const* a, unsigned n);

#endif
