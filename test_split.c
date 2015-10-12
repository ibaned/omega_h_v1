#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "element_field.h"
#include "global.h"
#include "inertia.h"
#include "ints.h"
#include "loop.h"
#include "files.h"
#include "mesh.h"
#include "subset.h"
#include "vtk.h"

static void split(struct mesh* m, char const* outpath, unsigned npieces,
    unsigned piece, unsigned depth)
{
  if (depth == 0) {
    if (piece == 0)
      write_pvtu(m, outpath, npieces, 0);
    char pvtupath[1024];
    char* filename;
    split_pathname(outpath, pvtupath, sizeof(pvtupath), &filename, 0);
    char prefix[1024];
    if (filename == pvtupath)
      strcpy(prefix, pvtupath);
    else
      sprintf(prefix, "%s/%s", pvtupath, filename);
    char vtupath[1024];
    parallel_pathname(prefix, npieces, piece, "vtu", vtupath, sizeof(vtupath));
    write_vtu(m, vtupath);
    free_mesh(m);
    return;
  }
  mesh_interp_to_elems(m, "coordinates");
  unsigned n = mesh_count(m, mesh_dim(m));
  double const* coords = mesh_find_tag(m, mesh_dim(m), "coordinates")->data;
  double* masses = 0;
  unsigned* in = local_inertia_mark(n, coords, masses);
  mesh_free_tag(m, mesh_dim(m), "coordinates");
  unsigned* offsets_in = uints_exscan(in, n);
  loop_free(in);
  struct mesh* m_in = subset_mesh(m, mesh_dim(m), offsets_in);
  split(m_in, outpath, npieces * 2, piece * 2 + 0, depth - 1);
  unsigned* offsets_out = uints_negate_offsets(offsets_in, n);
  loop_free(offsets_in);
  struct mesh* m_out = subset_mesh(m, mesh_dim(m), offsets_out);
  loop_free(offsets_out);
  split(m_out, outpath, npieces * 2, piece * 2 + 1, depth - 1);
  free_mesh(m);
}

int main(int argc, char** argv)
{
  assert(argc == 4);
  char const* inpath = argv[1];
  char const* outpath = argv[2];
  unsigned depth = (unsigned) atoi(argv[3]);
  struct mesh* m = read_vtu(inpath);
  mesh_number_simply(m);
  split(m, outpath, 1, 0, depth);
}
