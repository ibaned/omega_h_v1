#ifndef LOOP_H
#define LOOP_H

void* loop_malloc(unsigned long n);
void loop_free(void* p);

void* loop_host_malloc(unsigned long n);
void* loop_host_realloc(void* p, unsigned long n);
void loop_host_free(void* p);

#endif
