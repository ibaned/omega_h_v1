#ifndef LOOP_HOST_HPP
#define LOOP_HOST_HPP

namespace omega_h {

void* loop_host_malloc(unsigned long n);
#define LOOP_HOST_MALLOC(T, n) \
  static_cast<T*>(loop_host_malloc(sizeof(T) * (n)))
void* loop_host_realloc(void* p, unsigned long n);
#define LOOP_HOST_REALLOC(T, p, n) \
  static_cast<T*>(loop_host_realloc(p, sizeof(T) * (n)))
void loop_host_free(void* p);

unsigned loop_host_atomic_increment(unsigned* p);

void* loop_host_copy(void const* p, unsigned long n);
#define LOOP_HOST_COPY(T, p, n) \
  static_cast<T*>(loop_host_copy(p, sizeof(T) * (n)))

unsigned long loop_host_memory(void);
unsigned long loop_host_high_water(void);

}

#endif
