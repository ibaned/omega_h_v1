CC = gcc
CFLAGS = -std=c99 -g -O2 -fopenmp
LDFLAGS = -fno-omit-frame-pointer -fopenmp
LDLIBS = -lm
LOOP_MODE = openmp
