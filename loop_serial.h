#ifndef LOOP_SERIAL_H
#define LOOP_SERIAL_H

#include "loop_host.h"

#define LOOP_MALLOC(T, n) LOOP_HOST_MALLOC(T, n)
#define loop_free loop_host_free

#define LOOP_KERNEL(fname, ...) \
static void fname(__VA_ARGS__, unsigned i) \
{

#define LOOP_EXEC(fname, n, ...) \
for (unsigned i = 0; i < n; ++i) \
  fname(__VA_ARGS__, i);

#endif
