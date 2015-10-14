#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "base64.h"
#include "loop.h"

int main()
{
  char const* in = "any carnal pleasur";
  char* out = base64_encode(in, strlen(in));
  printf("%s\n", out);
  unsigned long size;
  void* outin = base64_decode(out, &size);
  assert(size == strlen(in));
  loop_host_free(out);
  fwrite(outin, 1, size, stdout);
  printf("\n");
  loop_host_free(outin);
}
