CXX = nvcc
CPPFLAGS = -x cu
CXXFLAGS = -g -O2 --gpu-code=sm_37 --gpu-architecture=compute_37 -rdc true
LDFLAGS = --gpu-code=sm_37 --gpu-architecture=compute_37 -rdc true
LOOP_MODE = cuda
