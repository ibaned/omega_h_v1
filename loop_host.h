#ifndef LOOP_HOST_H
#define LOOP_HOST_H

void* loop_host_malloc(unsigned long n);
void* loop_host_realloc(void* p, unsigned long n);
void loop_host_free(void* p);

#endif
