include /home/daibane/kokkos-install/Makefile.kokkos
CXX = nvcc
CXXFLAGS = -g -O2 -rdc true -expt-extended-lambda
CPPFLAGS = $(KOKKOS_CPPFLAGS) -x cu --gpu-code=sm_37 \
  --gpu-architecture=compute_37 --std=c++11
LDFLAGS = --gpu-code=sm_37 --gpu-architecture=compute_37 \
  $(KOKKOS_LDFLAGS)
LDLIBS = $(KOKKOS_LIBS)
LOOP_MODE = kokkos
