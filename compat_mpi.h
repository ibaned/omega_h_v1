#ifndef MPI_COMPAT_H
#define MPI_COMPAT_H

#include <mpi.h>

#include <assert.h>

#define CALL(f) do { int err = (f); assert(err == MPI_SUCCESS); } while(0)

struct comm {
  MPI_Comm c;
};

int compat_Neighbor_alltoallv(
    const void *sendbuf,
    const int sendcounts[],
    const int sdispls[],
    MPI_Datatype sendtype,
    void *recvbuf,
    const int recvcounts[],
    const int rdispls[],
    MPI_Datatype recvtype,
    MPI_Comm comm);

#endif
