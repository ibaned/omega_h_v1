#include "global.h"

#include "comm.h"
#include "ints.h"
#include "loop.h"

unsigned long* globalize_offsets(unsigned* local, unsigned n)
{
  unsigned long* global = LOOP_MALLOC(unsigned long, n);
  unsigned long lsum = local[n];
  unsigned long goff = comm_exscan_ulong(lsum);
  for (unsigned i = 0; i < n; ++i)
    global[i] = local[i] + goff;
  return global;
}

void global_to_linpart(unsigned long const* global, unsigned n,
    unsigned long total, unsigned nparts,
    unsigned** p_part, unsigned** p_local)
{
  unsigned long quot = total / nparts;
  unsigned long rem = total % nparts;
  unsigned* part = LOOP_MALLOC(unsigned, n);
  unsigned* local = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i) {
    unsigned long g = global[i];
    if (g < ((quot + 1) * rem)) {
      part[i] = (unsigned) (g / (quot + 1));
      local[i] = (unsigned) (g % (quot + 1));
    } else {
      g -= (quot + 1) * rem;
      part[i] = (unsigned) (g / quot + rem);
      local[i] = (unsigned) (g % quot);
    }
  }
  *p_part = part;
  *p_local = local;
}

unsigned linpart_size(unsigned long total, unsigned nparts, unsigned part)
{
  unsigned long quot = total / nparts;
  unsigned long rem = total % nparts;
  if (part < rem)
    return (unsigned) (quot + 1);
  return (unsigned) quot;
}
