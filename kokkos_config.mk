include /Users/dibanez/code/kokkos-install/Makefile.kokkos
CXX = clang++
CXXFLAGS = -g -O2 -Weverything -Werror -Wno-padded
CXXFLAGS += -Wno-c++98-compat -Wno-float-equal
CPPFLAGS = $(KOKKOS_CPPFLAGS) --std=c++11
LDFLAGS = $(KOKKOS_LDFLAGS)
LDLIBS = $(KOKKOS_LIBS)
LOOP_MODE = kokkos
