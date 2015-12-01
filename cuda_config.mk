CPP = nvcc
CC = nvcc
CFLAGS = -g -O2 -x cu --gpu-code=sm_37 --gpu-architecture=compute_37
LOOP_MODE = cuda
