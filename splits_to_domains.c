#include "splits_to_domains.h"

#include <assert.h>

#include "arrays.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

LOOP_KERNEL(split_to_domain,
    unsigned const* srcs_of_doms,
    unsigned srcs_per_dom,
    unsigned const* offset_of_srcs,
    unsigned const* vert_of_srcs,
    unsigned* dom_will_split,
    unsigned* direction_of_doms,
    unsigned* vert_of_doms)
  unsigned const* srcs_of_dom = srcs_of_doms + i * srcs_per_dom;
  for (unsigned j = 0; j < srcs_per_dom; ++j) {
    unsigned src = srcs_of_dom[j];
    if (offset_of_srcs[src] == offset_of_srcs[src + 1])
      continue;
    dom_will_split[i] = 1;
    direction_of_doms[i] = j;
    vert_of_doms[i] = vert_of_srcs[src];
    break;
  }
}

void project_splits_to_domains(
    unsigned dom_dim,
    unsigned src_dim,
    unsigned ndoms,
    unsigned const* srcs_of_doms,
    unsigned const* offset_of_srcs,
    unsigned const* vert_of_srcs,
    unsigned** p_offset_of_doms,
    unsigned** p_direction_of_doms,
    unsigned** p_vert_of_doms)
{
  assert(dom_dim >= src_dim);
  unsigned srcs_per_dom = the_down_degrees[dom_dim][src_dim];
  unsigned* dom_will_split = uints_filled(ndoms, 0);
  unsigned* direction_of_doms = LOOP_MALLOC(unsigned, ndoms);
  unsigned* vert_of_doms = LOOP_MALLOC(unsigned, ndoms);
  LOOP_EXEC(split_to_domain, ndoms,
      srcs_of_doms, srcs_per_dom, offset_of_srcs, vert_of_srcs,
      dom_will_split, direction_of_doms, vert_of_doms);
  *p_offset_of_doms = uints_exscan(dom_will_split, ndoms);
  loop_free(dom_will_split);
  *p_direction_of_doms = direction_of_doms;
  *p_vert_of_doms = vert_of_doms;
}

void mesh_splits_to_domains(
  struct mesh* m,
  unsigned dom_dim,
  unsigned src_dim,
  unsigned const* offset_of_srcs,
  unsigned const* vert_of_srcs,
  unsigned** p_offset_of_doms,
  unsigned** p_direction_of_doms,
  unsigned** p_vert_of_doms)
{
  unsigned ndoms = mesh_count(m, dom_dim);
  unsigned const* srcs_of_doms = mesh_ask_down(m, dom_dim, src_dim);
  project_splits_to_domains(dom_dim, src_dim, ndoms,
      srcs_of_doms, offset_of_srcs, vert_of_srcs,
      p_offset_of_doms, p_direction_of_doms, p_vert_of_doms);
}
