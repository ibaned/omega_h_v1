#ifndef FILES_H
#define FILES_H

void split_filename(char const* filename, char* buf,
    unsigned buf_size, char** prefix);
void parallel_filename(char const* prefix, unsigned npieces,
    unsigned piece, char const* suffix, char* buf, unsigned buf_size);

#endif
