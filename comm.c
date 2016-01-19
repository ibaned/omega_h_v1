#include "comm.h"

#include <assert.h>

#include "loop.h"

#if USE_MPI

#include "compat_mpi.h"

static struct comm world = { MPI_COMM_WORLD };
static struct comm self = { MPI_COMM_SELF };
static struct comm* using = &world;

void comm_init(void)
{
  CALL(MPI_Init(0,0));
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

struct comm* comm_graph(struct comm* c,
    unsigned nout, unsigned const* out, unsigned const* outweights)
{
  struct comm* c2 = LOOP_HOST_MALLOC(struct comm, 1);
  int n = (int) nout;
  int sources[1];
  MPI_Comm_rank(c->c, sources);
  int degrees[1] = {n};
  int* destinations = LOOP_HOST_MALLOC(int, nout);
  int* weights = LOOP_HOST_MALLOC(int, nout);
  for (unsigned i = 0; i < nout; ++i) {
    destinations[i] = (int) out[i];
    weights[i] = (int) outweights[i];
  }
  CALL(MPI_Dist_graph_create(c->c, 1, sources, degrees, destinations,
        weights, MPI_INFO_NULL, 0, &c2->c));
  loop_host_free(destinations);
  loop_host_free(weights);
  return c2;
}

struct comm* comm_graph_exact(struct comm* c,
    unsigned nin, unsigned const* in, unsigned const* inweights,
    unsigned nout, unsigned const* out, unsigned const* outweights)
{
  struct comm* c2 = LOOP_HOST_MALLOC(struct comm, 1);
  int indegree = (int) nin;
  int* sources = LOOP_HOST_MALLOC(int, nin);
  int* sourceweights = LOOP_HOST_MALLOC(int, nin);
  for (unsigned i = 0; i < nin; ++i) {
    sources[i] = (int) in[i];
    sourceweights[i] = (int) inweights[i];
  }
  int outdegree = (int) nout;
  int* destinations = LOOP_HOST_MALLOC(int, nout);
  int* destweights = LOOP_HOST_MALLOC(int, nout);
  for (unsigned i = 0; i < nout; ++i) {
    destinations[i] = (int) out[i];
    destweights[i] = (int) outweights[i];
  }
  CALL(MPI_Dist_graph_create_adjacent(c->c,
        indegree, sources, sourceweights,
        outdegree, destinations, destweights,
        MPI_INFO_NULL, 0, &c2->c));
  loop_host_free(sources);
  loop_host_free(sourceweights);
  loop_host_free(destinations);
  loop_host_free(destweights);
  return c2;

}

void comm_recvs(struct comm* c,
    unsigned* nin, unsigned** in, unsigned** inweights)
{
  int indegree, outdegree, weighted;
  CALL(MPI_Dist_graph_neighbors_count(c->c, &indegree, &outdegree, &weighted));
  assert(weighted != 0);
  *nin = (unsigned) indegree;
  unsigned nout = (unsigned) outdegree;
  int* sources = LOOP_HOST_MALLOC(int, *nin);
  int* sourceweights = LOOP_HOST_MALLOC(int, *nin);
  int* destinations = LOOP_HOST_MALLOC(int, nout);
  int* destweights = LOOP_HOST_MALLOC(int, nout);
  CALL(MPI_Dist_graph_neighbors(c->c, indegree, sources, sourceweights,
        outdegree, destinations, destweights));
  *in = LOOP_HOST_MALLOC(unsigned, *nin);
  *inweights = LOOP_HOST_MALLOC(unsigned, *nin);
  for (unsigned i = 0; i < *nin; ++i) {
    (*in)[i] = (unsigned) sources[i];
    (*inweights)[i] = (unsigned) sourceweights[i];
  }
  loop_host_free(sources);
  loop_host_free(sourceweights);
  loop_host_free(destinations);
  loop_host_free(destweights);
}

static void comm_exch_any(struct comm* c,
    unsigned width,
    void const* out, unsigned const* outcounts, unsigned const* outoffsets,
    void* in, unsigned const* incounts, unsigned const* inoffsets,
    MPI_Datatype type)
{
  int indegree, outdegree, weighted;
  CALL(MPI_Dist_graph_neighbors_count(c->c, &indegree, &outdegree, &weighted));
  int* sendcounts = LOOP_HOST_MALLOC(int, (unsigned) outdegree);
  int* sdispls = LOOP_HOST_MALLOC(int, (unsigned) outdegree);
  for (int i = 0; i < outdegree; ++i) {
    sendcounts[i] = (int) (outcounts[i] * width);
    sdispls[i] = (int) (outoffsets[i] * width);
  }
  int* recvcounts = LOOP_HOST_MALLOC(int, (unsigned) indegree);
  int* rdispls = LOOP_HOST_MALLOC(int, (unsigned) indegree);
  for (int i = 0; i < indegree; ++i) {
    recvcounts[i] = (int) (incounts[i] * width);
    rdispls[i] = (int) (inoffsets[i] * width);
  }
  CALL(compat_Neighbor_alltoallv(out, sendcounts, sdispls, type,
        in, recvcounts, rdispls, type, c->c));
  loop_host_free(sendcounts);
  loop_host_free(sdispls);
  loop_host_free(recvcounts);
  loop_host_free(rdispls);
}

void comm_exch_uints(struct comm* c,
    unsigned width,
    unsigned const* out, unsigned const* outcounts, unsigned const* outoffsets,
    unsigned* in, unsigned const* incounts, unsigned const* inoffsets)
{
  comm_exch_any(c, width, out, outcounts, outoffsets, in, incounts, inoffsets,
      MPI_UNSIGNED);
}

void comm_exch_doubles(struct comm* c,
    unsigned width,
    double const* out, unsigned const* outcounts, unsigned const* outoffsets,
    double* in, unsigned const* incounts, unsigned const* inoffsets)
{
  comm_exch_any(c, width, out, outcounts, outoffsets, in, incounts, inoffsets,
      MPI_DOUBLE);
}

void comm_exch_ulongs(struct comm* c,
    unsigned width,
    unsigned long const* out, unsigned const* outcounts, unsigned const* outoffsets,
    unsigned long* in, unsigned const* incounts, unsigned const* inoffsets)
{
  comm_exch_any(c, width, out, outcounts, outoffsets, in, incounts, inoffsets,
      MPI_UNSIGNED_LONG);
}

static void comm_sync_any(struct comm* c, void const* out, void* in, MPI_Datatype type)
{
  CALL(compat_Neighbor_allgather(out, 1, type, in, 1, type, c->c));
}

void comm_sync_uint(struct comm* c, unsigned out, unsigned* in)
{
  comm_sync_any(c, &out, in, MPI_UNSIGNED);
}

unsigned comm_bcast_uint(unsigned x)
{
  CALL(MPI_Bcast(&x, 1, MPI_UNSIGNED, 0, comm_using()->c));
  return x;
}

void comm_bcast_chars(char* s, unsigned n)
{
  CALL(MPI_Bcast(s, (int) n, MPI_CHAR, 0, comm_using()->c));
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
  CALL(MPI_Allreduce(MPI_IN_PLACE, p, (int) n, MPI_DOUBLE, MPI_SUM, using->c));
}

double comm_max_double(double x)
{
  CALL(MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_DOUBLE, MPI_MAX, using->c));
  return x;
}

unsigned long comm_add_ulong(unsigned long x)
{
  CALL(MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_UNSIGNED_LONG, MPI_SUM,
        using->c));
  return x;
}

unsigned long comm_exscan_ulong(unsigned long x)
{
  CALL(MPI_Exscan(MPI_IN_PLACE, &x, 1, MPI_UNSIGNED_LONG, MPI_SUM,
        using->c));
  if (!comm_rank())
    x = 0;
  return x;
}

unsigned long comm_max_ulong(unsigned long x)
{
  CALL(MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_UNSIGNED_LONG, MPI_MAX,
        using->c));
  return x;
}

unsigned comm_max_uint(unsigned x)
{
  CALL(MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_UNSIGNED, MPI_MAX,
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

struct split_comm {
  int dummy__;
};

struct comm* comm_split(struct comm* c, unsigned group, unsigned rank)
{
  (void)c;
  assert(group == 0);
  assert(rank == 0);
  return (struct comm*) LOOP_HOST_MALLOC(struct split_comm, 1);
}

struct graph_comm {
  unsigned nout;
  unsigned outweight;
};

struct comm* comm_graph(struct comm* c,
    unsigned nout, unsigned const* out, unsigned const* outweights)
{
  (void)c;
  struct graph_comm* gc = LOOP_HOST_MALLOC(struct graph_comm, 1);
  assert(nout <= 1);
  if (nout == 1)
    assert(out[0] == 0);
  gc->nout = nout;
  gc->outweight = outweights[0];
  return (struct comm*) gc;
}

struct comm* comm_graph_exact(struct comm* c,
    unsigned nin, unsigned const* in, unsigned const* inweights,
    unsigned nout, unsigned const* out, unsigned const* outweights)
{
  (void)c;
  struct graph_comm* gc = LOOP_HOST_MALLOC(struct graph_comm, 1);
  assert(nout <= 1);
  assert(nin == nout);
  if (nout == 1) {
    assert(out[0] == 0);
    assert(in[0] == 0);
    assert(inweights[0] == outweights[0]);
  }
  gc->nout = nout;
  gc->outweight = outweights[0];
  return (struct comm*) gc;
}

void comm_recvs(struct comm* c,
    unsigned* nin, unsigned** in, unsigned** incounts)
{
  struct graph_comm* gc = (struct graph_comm*) c;
  if (gc->nout == 0) {
    *nin = 0;
    *in = *incounts = 0;
  } else {
    assert(gc->nout == 1);
    *nin = 1;
    *in = LOOP_HOST_MALLOC(unsigned, 1);
    *incounts = LOOP_HOST_MALLOC(unsigned, 1);
    (*in)[0] = 0;
    /* assuming that the call to comm_graph was made
       with outweights = outcounts = incounts */
    (*incounts)[0] = gc->outweight;
  }
}

#define GENERIC_EXCH \
  (void) outoffsets; \
  (void) inoffsets; \
  struct graph_comm* gc = (struct graph_comm*) c; \
  if (gc->nout == 1) { \
    assert(outcounts[0] == incounts[0]); \
    for (unsigned i = 0; i < outcounts[0] * width; ++i) \
      in[i] = out[i]; \
  } \
  else \
    assert(gc->nout == 0);

void comm_exch_uints(struct comm* c,
    unsigned width,
    unsigned const* out, unsigned const* outcounts, unsigned const* outoffsets,
    unsigned* in, unsigned const* incounts, unsigned const* inoffsets)
{
  GENERIC_EXCH
}

void comm_exch_doubles(struct comm* c,
    unsigned width,
    double const* out, unsigned const* outcounts, unsigned const* outoffsets,
    double* in, unsigned const* incounts, unsigned const* inoffsets)
{
  GENERIC_EXCH
}

void comm_exch_ulongs(struct comm* c,
    unsigned width,
    unsigned long const* out, unsigned const* outcounts, unsigned const* outoffsets,
    unsigned long* in, unsigned const* incounts, unsigned const* inoffsets)
{
  GENERIC_EXCH
}

void comm_sync_uint(struct comm* c, unsigned out, unsigned* in)
{
  struct graph_comm* gc = (struct graph_comm*) c;
  if (gc->nout)
    in[0] = out;
}

unsigned comm_bcast_uint(unsigned x)
{
  return x;
}

void comm_bcast_chars(char* s, unsigned n)
{
  (void) s;
  (void) n;
}

void comm_free(struct comm* c)
{
  loop_host_free(c);
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

double comm_max_double(double x)
{
  return x;
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

unsigned comm_max_uint(unsigned x)
{
  return x;
}

#endif

double comm_add_double(double x)
{
  double a[1] = {x};
  comm_add_doubles(a, 1);
  return a[0];
}
