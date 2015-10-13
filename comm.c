#include "comm.h"

#if USE_MPI

#include <assert.h>
#include <mpi.h>

#include "loop.h"

struct comm {
  MPI_Comm c;
};

#define CALL(f) do { int err = (f); assert(err == MPI_SUCCESS); } while(0)

static struct comm world = { MPI_COMM_WORLD };
static struct comm self = { MPI_COMM_SELF };
static struct comm* using = &world;

void comm_init(void)
{
  CALL(MPI_Init(NULL,NULL));
}

void comm_fini(void)
{
  CALL(MPI_Finalize());
}

struct comm* comm_world(void)
{
  return &world;
}

struct comm* comm_self(void)
{
  return &self;
}

struct comm* comm_using(void)
{
  return using;
}

void comm_use(struct comm* c)
{
  using = c;
}

struct comm* comm_split(struct comm* c, unsigned group, unsigned rank)
{
  struct comm* c2 = LOOP_HOST_MALLOC(struct comm, 1);
  CALL(MPI_Comm_split(c->c, (int) group, (int) rank, &c2->c));
  return c2;
}

void comm_free(struct comm* c)
{
  assert(c != &world);
  assert(c != &self);
  CALL(MPI_Comm_free(&c->c));
  loop_host_free(c);
}

unsigned comm_rank(void)
{
  int rank;
  CALL(MPI_Comm_rank(using->c, &rank));
  return (unsigned) rank;
}

unsigned comm_size(void)
{
  int size;
  CALL(MPI_Comm_size(using->c, &size));
  return (unsigned) size;
}

void comm_add_doubles(double* p, unsigned n)
{
  CALL(MPI_Allreduce(p, MPI_IN_PLACE, (int) n, MPI_DOUBLE, MPI_SUM, using->c));
}

unsigned long comm_add_ulong(unsigned long x)
{
  CALL(MPI_Allreduce(&x, MPI_IN_PLACE, 1, MPI_UNSIGNED_LONG, MPI_SUM,
        using->c));
  return x;
}

unsigned long comm_exscan_ulong(unsigned long x)
{
  CALL(MPI_Exscan(&x, MPI_IN_PLACE, 1, MPI_UNSIGNED_LONG, MPI_SUM,
        using->c));
  if (!comm_rank())
    x = 0;
  return x;
}

unsigned long comm_max_ulong(unsigned long x)
{
  CALL(MPI_Allreduce(&x, MPI_IN_PLACE, 1, MPI_UNSIGNED_LONG, MPI_MAX,
        using->c));
  return x;
}

#else

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

unsigned long comm_max_ulong(unsigned long x)
{
  return x;
}

#endif
