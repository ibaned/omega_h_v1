#ifndef LOOP_SERIAL_H
#define LOOP_SERIAL_H

#include "loop_host.h"

#define loop_malloc loop_host_malloc
#define loop_free loop_host_free

#define LOOP_KERNEL(fname, ...) \
static void fname(__VA_ARGS__, unsigned i) \
{

#define LOOP_EXEC(fname, n, ...) \
for (unsigned i = 0; i < n; ++i) \
  fname(__VA_ARGS__, i);

#endif
