include /home1/03049/ibaned/kokkos-install/Makefile.kokkos
CXX = mpicxx
CXXFLAGS = -g -O2 -fopenmp -mmic
CPPFLAGS = $(KOKKOS_CPPFLAGS) --std=c++11
LDFLAGS = $(KOKKOS_LDFLAGS)
LDLIBS = $(KOKKOS_LIBS)
LOOP_MODE = kokkos
