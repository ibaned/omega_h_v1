include /home/daibane/kokkos-install/Makefile.kokkos
CXX = clang++
CXXFLAGS = -g -O2 -Weverythig -Werror -Wno-padded
CPPFLAGS = $(KOKKOS_CPPFLAGS) --std=c++11
LDFLAGS = $(KOKKOS_LDFLAGS)
LDLIBS = $(KOKKOS_LIBS)
LOOP_MODE = kokkos
