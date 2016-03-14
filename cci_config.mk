#module load gnu-4.7.2
CXX = mpicxx
CXXFLAGS = -O2 -Wall -Wextra -pedantic -std=c++11
LDLIBS = -lm
USE_MPI = 1
# CCI's MPI library is older than MPI3.
USE_MPI3 = 0
