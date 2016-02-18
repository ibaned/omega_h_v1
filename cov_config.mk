CC = mpicc
CFLAGS = -g -O0 -std=c99 -fprofile-arcs -ftest-coverage
USE_ZLIB = 1
USE_MPI = 1
USE_MPI3 = 1
PATIENT = 1
LDFLAGS = -fprofile-arcs
LDLIBS = -lz
