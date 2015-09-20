#include "comm.h"

void comm_init(void)
{
}

void comm_fini(void)
{
}

struct comm* comm_world(void)
{
  return (struct comm*)1;
}

struct comm* comm_self(void)
{
  return (struct comm*)1;
}

unsigned comm_rank(struct comm* c)
{
  (void)c;
  return 0;
}

unsigned comm_size(struct comm* c)
{
  (void)c;
  return 1;
}

void comm_add_doubles(struct comm* c, double* p, unsigned n)
{
  (void) c;
  (void) p;
  (void) n;
}

unsigned long comm_add_ulong(struct comm* c, unsigned long x)
{
  (void) c;
  return x;
}
