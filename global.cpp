#include "global.hpp"

#include "arrays.hpp"
#include "comm.hpp"
#include "ints.hpp"
#include "int_casts.hpp"
#include "loop.hpp"

namespace omega_h {

LOOP_KERNEL(globalize_offset,
    unsigned const* local,
    unsigned long offset,
    unsigned long* global)
  global[i] = local[i] + offset;
}

unsigned long* globalize_offsets(unsigned* local, unsigned n)
{
  unsigned long* global = LOOP_MALLOC(unsigned long, n);
  unsigned long lsum = array_at(local, n);
  unsigned long offset = comm_exscan_ulong(lsum);
  LOOP_EXEC(globalize_offset, n, local, offset, global);
  return global;
}

LOOP_KERNEL(global_to_lin,
    unsigned long const* global,
    unsigned long quot,
    unsigned long rem,
    unsigned* local,
    unsigned* part)
  unsigned long g = global[i];
  if (g < ((quot + 1) * rem)) {
    part[i] = U(g / (quot + 1));
    local[i] = U(g % (quot + 1));
  } else {
    g -= (quot + 1) * rem;
    part[i] = U(g / quot + rem);
    local[i] = U(g % quot);
  }
}

void global_to_linpart(unsigned long const* global, unsigned n,
    unsigned long total, unsigned nparts,
    unsigned** p_part, unsigned** p_local)
{
  unsigned long quot = total / nparts;
  unsigned long rem = total % nparts;
  unsigned* part = LOOP_MALLOC(unsigned, n);
  unsigned* local = LOOP_MALLOC(unsigned, n);
  LOOP_EXEC(global_to_lin, n, global, quot, rem,
      local, part);
  *p_part = part;
  *p_local = local;
}

unsigned linpart_size(unsigned long total, unsigned nparts, unsigned part)
{
  unsigned long quot = total / nparts;
  unsigned long rem = total % nparts;
  if (part < rem)
    return U(quot + 1);
  return U(quot);
}

}
