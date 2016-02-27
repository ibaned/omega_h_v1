CPP = nvcc
CC = nvcc
CPPFLAGS = -x cu
CFLAGS = -g -O2 --gpu-code=sm_30 --gpu-architecture=compute_30 -rdc true
LDFLAGS = --gpu-code=sm_30 --gpu-architecture=compute_30 -rdc true
LOOP_MODE = cuda
