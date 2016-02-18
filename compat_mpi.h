#ifndef MPI_COMPAT_H
#define MPI_COMPAT_H

#include <mpi.h>

#if MPI_VERSION < 2 || (MPI_VERSION == 2 && MPI_SUBVERSION < 2)
#error "omega_h does not support MPI versions less than 2.2"
#endif

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

int compat_Neighbor_allgather(
    const void *sendbuf,
    int sendcount,
    MPI_Datatype sendtype,
    void *recvbuf,
    int recvcount,
    MPI_Datatype recvtype,
    MPI_Comm comm);

#endif
