#include "sliver_keys.h"
#include "tables.h"
#include "quality.h"
#include "algebra.h"
#include <stdlib.h>

void sliver_keys(
    unsigned elem_dim,
    unsigned nelems,
    unsigned const* verts_of_elems,
    double const* coords,
    enum sliver_type target,
    double qual_floor,
    double edge_ratio_floor,
    unsigned** bad_elems_out,
    unsigned** key_of_elems_out)
{
  unsigned* bad_elems = malloc(sizeof(unsigned) * nelems);
  unsigned* key_of_elems = malloc(sizeof(unsigned) * nelems);
  sliver_type_function stf = the_sliver_type_functions[elem_dim];
  unsigned verts_per_elem = the_down_degrees[elem_dim][0];
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* verts_of_elem = verts_of_elems + i * verts_per_elem;
    double elem_coords[MAX_DOWN][3];
    for (unsigned j = 0; j < verts_per_elem; ++j) {
      unsigned vert = verts_of_elem[j];
      copy_vector(coords + vert * 3, elem_coords[j], 3);
    }
    unsigned key;
    enum sliver_type st = stf(elem_coords, qual_floor, edge_ratio_floor, &key);
    if (st != target)
      bad_elems[i] = 0;
    else {
      bad_elems[i] = 1;
      key_of_elems[i] = key;
    }
  }
  *bad_elems_out = bad_elems;
  *key_of_elems_out = key_of_elems;
}
