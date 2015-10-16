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
  char const* tmp = out;
  void* outin = base64_decode(&tmp, strlen(in));
  loop_host_free(out);
  fwrite(outin, 1, strlen(in), stdout);
  printf("\n");
  loop_host_free(outin);
}
