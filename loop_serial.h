#ifndef LOOP_SERIAL_H
#define LOOP_SERIAL_H

#include "loop_host.h"

#define LOOP_MALLOC(T, n) LOOP_HOST_MALLOC(T, n)
#define loop_free loop_host_free

#define loop_to_host loop_host_copy
#define loop_to_device loop_host_copy
#define loop_memcpy loop_host_memcpy

#define loop_atomic_increment loop_host_atomic_increment

#define LOOP_KERNEL(fname, ...) \
static void fname(__VA_ARGS__, unsigned i) \
{

#define LOOP_EXEC(fname, n, ...) \
for (unsigned i = 0; i < n; ++i) \
  fname(__VA_ARGS__, i);

unsigned loop_size(void);

#endif
