#include "compat_mpi.h"

#include "loop_host.h"

/* this file exists to support MPI_VERSION < 3.
   it implements some of the high-level Neighbor
   functions used by omega_h */

#if MPI_VERSION >= 3 && USE_MPI3

int compat_Neighbor_alltoallv(
    const void *sendbuf,
    const int sendcounts[],
    const int sdispls[],
    MPI_Datatype sendtype,
    void *recvbuf,
    const int recvcounts[],
    const int rdispls[],
    MPI_Datatype recvtype,
    MPI_Comm comm)
{
  return MPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype,
      recvbuf, recvcounts, rdispls, recvtype, comm);
}

int compat_Neighbor_allgather(
    const void *sendbuf,
    int sendcount,
    MPI_Datatype sendtype,
    void *recvbuf,
    int recvcount,
    MPI_Datatype recvtype,
    MPI_Comm comm)
{
  return MPI_Neighbor_allgather(sendbuf, sendcount, sendtype,
      recvbuf, recvcount, recvtype, comm);
}

#else

#define MY_TAG 42

int compat_Neighbor_alltoallv(
    const void *sendbuf,
    const int sendcounts[],
    const int sdispls[],
    MPI_Datatype sendtype,
    void *recvbuf,
    const int recvcounts[],
    const int rdispls[],
    MPI_Datatype recvtype,
    MPI_Comm comm)
{
  int indegree, outdegree, weighted;
  CALL(MPI_Dist_graph_neighbors_count(comm, &indegree, &outdegree, &weighted));
  int* sources = LOOP_HOST_MALLOC(int, (unsigned) indegree);
  int* sourceweights = LOOP_HOST_MALLOC(int, (unsigned) indegree);
  int* destinations = LOOP_HOST_MALLOC(int, (unsigned) outdegree);
  int* destweights = LOOP_HOST_MALLOC(int, (unsigned) outdegree);
  CALL(MPI_Dist_graph_neighbors(comm, indegree, sources, sourceweights,
        outdegree, destinations, destweights));
  loop_host_free(sourceweights);
  loop_host_free(destweights);
  int sendwidth;
  CALL(MPI_Type_size(sendtype, &sendwidth));
  int recvwidth;
  CALL(MPI_Type_size(sendtype, &recvwidth));
  MPI_Request* recvreqs = LOOP_HOST_MALLOC(MPI_Request, (unsigned) indegree);
  MPI_Request* sendreqs = LOOP_HOST_MALLOC(MPI_Request, (unsigned) outdegree);
  for (int i = 0; i < indegree; ++i)
    CALL(MPI_Irecv(((char*)recvbuf) + rdispls[i] * recvwidth,
          recvcounts[i], recvtype, sources[i], MY_TAG, comm,
          recvreqs + i));
  loop_host_free(sources);
  /* putting a barrier here may increase performance by ensuring
     that receives are posted before their corresponding sends,
     but of course decreases performance because there is a barrier.
     add/remove this line as you wish. */
  CALL(MPI_Barrier(comm));
  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(((char const*)sendbuf) + sdispls[i] * sendwidth,
          sendcounts[i], sendtype, destinations[i], MY_TAG, comm,
          sendreqs + i));
  loop_host_free(destinations);
  CALL(MPI_Waitall(outdegree, sendreqs, MPI_STATUSES_IGNORE));
  loop_host_free(sendreqs);
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  loop_host_free(recvreqs);
  return MPI_SUCCESS;
}

int compat_Neighbor_allgather(
    const void *sendbuf,
    int sendcount,
    MPI_Datatype sendtype,
    void *recvbuf,
    int recvcount,
    MPI_Datatype recvtype,
    MPI_Comm comm)
{
  int indegree, outdegree, weighted;
  CALL(MPI_Dist_graph_neighbors_count(comm, &indegree, &outdegree, &weighted));
  int* sources = LOOP_HOST_MALLOC(int, (unsigned) indegree);
  int* sourceweights = LOOP_HOST_MALLOC(int, (unsigned) indegree);
  int* destinations = LOOP_HOST_MALLOC(int, (unsigned) outdegree);
  int* destweights = LOOP_HOST_MALLOC(int, (unsigned) outdegree);
  CALL(MPI_Dist_graph_neighbors(comm, indegree, sources, sourceweights,
        outdegree, destinations, destweights));
  loop_host_free(sourceweights);
  loop_host_free(destweights);
  int recvwidth;
  CALL(MPI_Type_size(sendtype, &recvwidth));
  MPI_Request* recvreqs = LOOP_HOST_MALLOC(MPI_Request, (unsigned) indegree);
  MPI_Request* sendreqs = LOOP_HOST_MALLOC(MPI_Request, (unsigned) outdegree);
  for (int i = 0; i < indegree; ++i)
    CALL(MPI_Irecv(((char*)recvbuf) + i * recvcount * recvwidth,
          recvcount, recvtype, sources[i], MY_TAG, comm,
          recvreqs + i));
  loop_host_free(sources);
  CALL(MPI_Barrier(comm));
  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(sendbuf,
          sendcount, sendtype, destinations[i], MY_TAG, comm,
          sendreqs + i));
  loop_host_free(destinations);
  CALL(MPI_Waitall(outdegree, sendreqs, MPI_STATUSES_IGNORE));
  loop_host_free(sendreqs);
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  loop_host_free(recvreqs);
  return MPI_SUCCESS;
}

#endif
