#module load gnu-4.7.2
CC = mpicc
CPP = mpicc
CFLAGS = -O2 -Wall -Wextra -pedantic -std=c99
LDFLAGS = -O2
LDLIBS = -lm
USE_MPI = 1
# CCI's MPI library is older than MPI3.
USE_MPI3 = 0
