#include "comm.h"

struct comm;

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

struct comm* comm_using(void)
{
  return (struct comm*)1;
}

void comm_use(struct comm* c)
{
  (void)c;
}

struct comm* comm_split(struct comm* c, unsigned group, unsigned rank)
{
  (void)c;
  (void)group;
  (void)rank;
  return (struct comm*)1;
}

void comm_free(struct comm* c)
{
  (void)c;
}

unsigned comm_rank(void)
{
  return 0;
}

unsigned comm_size(void)
{
  return 1;
}

void comm_add_doubles(double* p, unsigned n)
{
  (void)p;
  (void)n;
}

unsigned long comm_add_ulong(unsigned long x)
{
  return x;
}

unsigned long comm_exscan_ulong(unsigned long x)
{
  (void)x;
  return 0;
}
