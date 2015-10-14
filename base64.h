#ifndef BASE_64_H
#define BASE_64_H

char* base64_encode(void const* data, unsigned long size);
void* base64_decode(char const* text, unsigned long* size);

void print_base64_reverse(void);

#endif
