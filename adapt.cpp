#include "adapt.hpp"

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "coarsen.hpp"
#include "comm.hpp"
#include "doubles.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "quality.hpp"
#include "refine.hpp"
#include "size.hpp"
#include "swap.hpp"

namespace omega_h {

static unsigned global_op_count = 0;
static unsigned global_max_ops = 0;

void adapt_summary(struct mesh* m)
{
  unsigned long total_elems = comm_add_ulong(mesh_count(m, mesh_dim(m)));
  double minqual = comm_min_double(mesh_min_quality(m));
  unsigned nedges = mesh_count(m, 1);
  double* edge_sizes = mesh_measure_edges_for_adapt(m);
  double min = comm_min_double(doubles_min(edge_sizes, nedges));
  double max = comm_max_double(doubles_max(edge_sizes, nedges));
  loop_free(edge_sizes);
  if (comm_rank() == 0)
    printf("%10lu elements, min quality %.0f%%, metric range %.2f - %.2f\n",
        total_elems, minqual * 100.0, min, max);
}

static void incr_op_count(struct mesh* m)
{
  if (global_op_count > global_max_ops) {
    fprintf(stderr, "mesh_adapt could not succeed after %u operations\n",
        global_max_ops);
    abort();
  }
  ++global_op_count;
  adapt_summary(m);
}

static void satisfy_size(struct mesh* m, double size_floor, double good_qual)
{
  double qual_floor = mesh_min_quality(m);
  if (good_qual < qual_floor)
    qual_floor = good_qual;
  while (refine_by_size(m, qual_floor))
    incr_op_count(m);
  while (coarsen_by_size(m, qual_floor, size_floor)) {
    fprintf(stderr, "coarsen_by_size returned nonzero, calling incr_op_count\n");
    incr_op_count(m);
  }
}

static void satisfy_shape(
    struct mesh* m,
    double qual_floor,
    unsigned nsliver_layers)
{
  while (1) {
    double prev_qual = comm_min_double(mesh_min_quality(m));
    if (prev_qual >= qual_floor)
      return;
    if (mesh_dim(m) == 3 &&
        swap_slivers(m, qual_floor, nsliver_layers)) {
      incr_op_count(m);
      continue;
    }
    if (coarsen_slivers(m, qual_floor, nsliver_layers)) {
      fprintf(stderr, "coarsen_slivers returned nonzero, calling incr_op_count\n");
      incr_op_count(m);
      continue;
    }
    break; //ran out of options, but return from adapt anyway
  }
}

unsigned mesh_adapt(struct mesh* m,
    double size_ratio_floor,
    double good_qual,
    unsigned nsliver_layers,
    unsigned max_ops)
{
  global_op_count = 0;
  global_max_ops = max_ops;
  adapt_summary(m);
  satisfy_size(m, size_ratio_floor, good_qual);
  satisfy_shape(m, good_qual, nsliver_layers);
  return global_op_count > 0;
}

}
