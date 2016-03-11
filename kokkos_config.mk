CXX = mpicxx
CXXFLAGS = -g -O2 -Werror -Weverything -Wno-c++98-compat-pedantic -Wno-float-equal -Wno-padded
KOKKOS_PREFIX = /Users/dibanez/code/kokkos-install
LOOP_MODE = kokkos
USE_ZLIB = 1
USE_MPI = 1
USE_MPI3 = 1
PATIENT = 1
LDLIBS = -lz
VALGRIND = "valgrind --leak-check=full --error-exitcode=1"
CPPFLAGS += -std=c++11
CPPFLAGS += -I$(KOKKOS_PREFIX)/include
LDFLAGS += -L$(KOKKOS_PREFIX)/lib
LDLIBS += -lkokkos
