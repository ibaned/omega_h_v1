#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

namespace omega_h {

typedef std::chrono::time_point<std::chrono::high_resolution_clock> Now;

static inline Now now()
{
  return std::chrono::high_resolution_clock::now();
}

static inline double seconds_between(Now a, Now b)
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() * 1e-9;
}

}

#endif

