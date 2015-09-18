#ifndef COMM_H
#define COMM_H

struct comm;

void comm_init(void);
void comm_fini(void);
struct comm* comm_world(void);
unsigned comm_rank(struct comm* c);
unsigned comm_size(struct comm* c);
void comm_add_doubles(struct comm* c, double* p, unsigned n);
unsigned long comm_add_ulong(struct comm* c, unsigned long x);

#endif
