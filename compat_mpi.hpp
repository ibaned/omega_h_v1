#ifndef MPI_COMPAT_HPP
#define MPI_COMPAT_HPP

#ifdef __clang__
#pragma clang system_header
#endif

#ifdef __bgq__
/*
The BlueGene/Q MPI headers have MPICH2_CONST
everywhere that const should be, and apparently
define it to nothing by default, which causes
compile errors for code that actually is const-correct.
so, "#define MPICH2_CONST const" would be enough. but no.
They also have this line of code in mpicxx.h:
    virtual Datatype Create_indexed( int v1, const MPICH2_CONST int * v2,
                                      const MPICH2_CONST int * v3 ) const
Which triggers a compile error about duplicate const.
*/
#define MPICH2_CONST const
#define MPICH_SKIP_MPICXX
#endif

#include <mpi.h>

#if MPI_VERSION < 2 || (MPI_VERSION == 2 && MPI_SUBVERSION < 2)
#error "omega_h does not support MPI versions less than 2.2"
#endif

#include <cassert>

#define CALL(f) do { int err = (f); assert(err == MPI_SUCCESS); } while(0)

namespace omega_h {

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

}

#endif
