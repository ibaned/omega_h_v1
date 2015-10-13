#include "global.h"

#include "comm.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tag.h"

void mesh_number_simply(struct mesh* m)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned long* data = LOOP_MALLOC(unsigned long, nverts);
  for (unsigned i = 0; i < nverts; ++i)
    data[i] = i;
  mesh_add_tag(m, 0, TAG_U64, "global_number", 1, data);
}

unsigned long* globalize_offsets(unsigned* local, unsigned n)
{
  unsigned long* global = LOOP_MALLOC(unsigned long, n);
  unsigned long lsum = local[n];
  unsigned long goff = comm_exscan_ulong(lsum);
  for (unsigned i = 0; i < n; ++i)
    global[i] = local[i] + goff;
  return global;
}

void global_to_parts(unsigned long* global, unsigned n,
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
