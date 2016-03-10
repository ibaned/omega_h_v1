#ifndef LOOP_SERIAL_HPP
#define LOOP_SERIAL_HPP

#include <assert.h>

#include "loop_host.hpp"

#define LOOP_MALLOC(T, n) LOOP_HOST_MALLOC(T, n)
#define loop_free loop_host_free

#define loop_atomic_increment loop_host_atomic_increment

unsigned loop_size(void);

#define LOOP_IN
#define LOOP_INOUT

#define LOOP_CONST

#define LOOP_NORETURN(x) assert(0)

#endif
