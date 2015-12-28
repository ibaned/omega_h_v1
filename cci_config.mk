#module load gnu-4.7.2
CC = mpicc
CPP = mpicc
CFLAGS = -O2 -Wall -Wextra -pedantic -std=c99
LDFLAGS = -flto -O2
LDLIBS = -lm
# CCI's MPI library is older than MPI3.
USE_MPI = 0
