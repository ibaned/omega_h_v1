#ifndef FILES_H
#define FILES_H

void split_pathname(char const* pathname, char* buf,
    unsigned buf_size, char** filename, char** suffix);
void parallel_pathname(char const* prefix, unsigned npieces,
    unsigned piece, char const* suffix, char* buf, unsigned buf_size);

#endif
