CXX = mpicxx
CXXFLAGS = -g -O2 -std=c++11
USE_ZLIB = 1
USE_MPI = 1
USE_MPI3 = 1
PATIENT = 1
LDLIBS = -lz
VALGRIND = "valgrind --leak-check=full --error-exitcode=1"
