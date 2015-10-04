#ifndef COMM_H
#define COMM_H

struct comm;

void comm_init(void);
void comm_fini(void);

struct comm* comm_world(void);
struct comm* comm_self(void);

struct comm* comm_using(void);
void comm_use(struct comm* c);

struct comm* comm_split(struct comm* c, unsigned group, unsigned rank);
void comm_free(struct comm* c);

unsigned comm_rank(void);
unsigned comm_size(void);

void comm_add_doubles(double* p, unsigned n);
unsigned long comm_add_ulong(unsigned long x);
unsigned long comm_max_ulong(unsigned long x);

#endif
