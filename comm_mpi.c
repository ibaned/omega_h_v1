#include "comm.h"
#include "loop.h"
#include <mpi.h>
#include <assert.h>

struct comm {
  MPI_Comm c;
};

#define CALL(f) do { int err = (f); assert(err == MPI_SUCCESS); } while(0)

void comm_init(void)
{
  CALL(MPI_Init(NULL,NULL));
}

void comm_fini(void)
{
  CALL(MPI_Finalize());
}

static struct comm world = { MPI_COMM_WORLD };

struct comm* comm_world(void)
{
  return &world;
}

unsigned comm_rank(struct comm* c)
{
  int rank;
  CALL(MPI_Comm_rank(c->c, &rank));
  return (unsigned) rank;
}

unsigned comm_size(struct comm* c)
{
  int size;
  CALL(MPI_Comm_size(c->c, &size));
  return (unsigned) size;
}

void comm_add_doubles(struct comm* c, double* p, unsigned n)
{
  CALL(MPI_Allreduce(p, MPI_IN_PLACE, (int) n, MPI_DOUBLE, MPI_SUM, c->c));
}

unsigned long comm_add_ulong(struct comm* c, unsigned long x)
{
  CALL(MPI_Allreduce(&x, MPI_IN_PLACE, 1, MPI_UNSIGNED_LONG, MPI_SUM, c->c));
  return x;
}
