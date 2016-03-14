#ifndef LOOP_SERIAL_HPP
#define LOOP_SERIAL_HPP

#include <cassert>

#include "loop_host.hpp"

namespace omega_h {

#define LOOP_MALLOC(T, n) \
  static_cast<T*>(loop_host_malloc(sizeof(T) * (n)))
#define loop_free loop_host_free

#define loop_atomic_increment loop_host_atomic_increment

#define LOOP_KERNEL(fname, ...) \
static void fname(__VA_ARGS__, unsigned i) \
{

#define LOOP_EXEC(fname, n, ...) \
for (unsigned i = 0; i < n; ++i) \
  fname(__VA_ARGS__, i);

#define LOOP_IN
#define LOOP_INOUT

#define LOOP_CONST

#define LOOP_NORETURN(x) assert(0)

}

#endif
