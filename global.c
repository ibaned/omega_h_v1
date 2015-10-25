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

void global_to_linpart(unsigned long* global, unsigned n,
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

/* given an array that indicates which part an
   entry belongs to (the unique set of parts is
   an arbitrary but usually small set of numbers),
   this function assigns a contiguous ID
   to each part (i.e. numbering the set of unique
   entires in (parts)).
   it also numbers the entries belonging to each
   part consecutively.
   the resulting category (cats) and category
   index (cat_indices) arrays can be used
   so sort the input data by category in linear time
   (i.e. this function is a generalized sort).
   we rely on the assumption that the number of
   unique parts is small to settle on a runtime
   that is O(n * (# unique parts) * log(n)).
   the log(n) comes from the runtime of uints_exscan.
*/

void categorize_by_part(unsigned const* parts, unsigned n,
    unsigned** p_cats, unsigned** p_cat_indices,
    unsigned* p_ncats,
    unsigned** p_cat_parts, unsigned** p_cat_counts)
{
  /* queued[i]==1 iff (i) has not been categorized yet */
  unsigned* queued = LOOP_MALLOC(unsigned, n);
  for (unsigned i = 0; i < n; ++i)
    queued[i] = 1;
  /* loop over categories, we don't know how many but
     there certainly won't be more than (n) */
  unsigned* cats = LOOP_MALLOC(unsigned, n);
  unsigned* cat_indices = LOOP_MALLOC(unsigned, n);
  unsigned* cat_counts = LOOP_MALLOC(unsigned, n);
  unsigned* cat_parts = LOOP_MALLOC(unsigned, n);
  unsigned cat;
  for (cat = 0; cat < n; ++cat) {
    unsigned current_part = 0;
    unsigned* queue_offsets = uints_exscan(queued, n);
    unsigned nqueued = queue_offsets[n];
    if (nqueued == 0) {
      loop_free(queue_offsets);
      break; /* stop when all entries categorized */
    }
    /* process the part of the first uncategorized entry */
    for (unsigned i = 0; i < n; ++i)
      if ((queue_offsets[i + 1] - queue_offsets[i] == 1) &&
          queue_offsets[i] == 0)
        current_part = parts[i];
    cat_parts[cat] = current_part;
    loop_free(queue_offsets);
    unsigned* in_part = LOOP_MALLOC(unsigned, n);
    for (unsigned i = 0; i < n; ++i) {
      if (parts[i] == current_part) {
        cats[i] = cat;
        in_part[i] = 1;
        queued[i] = 0;
      } else {
        in_part[i] = 0;
      }
    }
    unsigned* part_indices = uints_exscan(in_part, n);
    cat_counts[cat] = part_indices[n];
    for (unsigned i = 0; i < n; ++i)
      if (in_part[i])
        cat_indices[i] = part_indices[i];
    loop_free(in_part);
    loop_free(part_indices);
  }
  unsigned ncats = cat;
  loop_free(queued);
  *p_cats = cats;
  *p_cat_indices = cat_indices;
  *p_ncats = ncats;
  /* shrink these arrays to fit, they were
     allocated to the maximum possible size of (n) */
  *p_cat_parts = uints_copy(cat_parts, ncats);
  loop_free(cat_parts);
  *p_cat_counts = uints_copy(cat_counts, ncats);
  loop_free(cat_counts);
}

#define GENERIC_SORT_BY_CATEGORY(T, a, n, cats, indices, cat_offsets) \
  T* out = LOOP_MALLOC(T, n); \
  for (unsigned i = 0; i < n; ++i) { \
    unsigned cat = cats[i]; \
    unsigned idx = indices[i]; \
    unsigned o = cat_offsets[cat]; \
    out[o + idx] = a[i]; \
  } \
  return out;

unsigned* sort_uints_by_category(
    unsigned const* a,
    unsigned n,
    unsigned const* cats,
    unsigned const* indices,
    unsigned const* cat_offsets)
{
  GENERIC_SORT_BY_CATEGORY(unsigned, a, n, cats, indices, cat_offsets);
}
