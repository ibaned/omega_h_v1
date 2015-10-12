#include "global.h"

#include "comm.h"
#include "loop.h"
#include "mesh.h"
#include "tag.h"

void mesh_number_simply(struct mesh* m)
{
  unsigned nverts = mesh_count(m, 0);
  unsigned long* data = loop_malloc(sizeof(unsigned long) * nverts);
  for (unsigned i = 0; i < nverts; ++i)
    data[i] = i;
  mesh_add_tag(m, 0, TAG_U64, "global_number", 1, data);
}

unsigned long* globalize_offsets(unsigned* local, unsigned n)
{
  unsigned long* global = loop_malloc(sizeof(unsigned long) * n);
  unsigned long lsum = local[n];
  unsigned long goff = comm_exscan_ulong(lsum);
  for (unsigned i = 0; i < n; ++i)
    global[i] = local[i] + goff;
  return global;
}
