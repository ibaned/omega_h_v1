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
