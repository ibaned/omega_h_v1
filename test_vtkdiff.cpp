#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "include/omega_h.h"
#include "mesh.hpp"
#include "mesh_diff.hpp"
#include "vtk_io.hpp"

int main(int argc, char** argv)
{
  osh_init(&argc, &argv);
  double tol = 1e-6;
  double floor = 0;
  unsigned get_tol = 0;
  unsigned get_floor = 0;
  unsigned get_help = 0;
  char const* filea = 0;
  char const* fileb = 0;
  unsigned allow_superset = 0;
  for (int i = 1; i < argc; ++i) {
    if (get_tol) {
      tol = atof(argv[i]);
      get_tol = 0;
      continue;
    }
    if (get_floor) {
      floor = atof(argv[i]);
      get_floor = 0;
      continue;
    }
    if (!strcmp(argv[i], "-tolerance")) {
      get_tol = 1;
      continue;
    }
    if (!strcmp(argv[i], "-Floor")) {
      get_floor = 1;
      continue;
    }
    if (!strcmp(argv[i], "-help")) {
      get_help = 1;
      continue;
    }
    if (!strcmp(argv[i], "-superset")) {
      allow_superset = 1;
      continue;
    }
    if (!filea) {
      filea = argv[i];
      continue;
    }
    if (!fileb) {
      fileb = argv[i];
      continue;
    }
  }
  if (get_tol) {
    printf("-tolerance needs an argument\n");
    return -1;
  }
  if (get_floor) {
    printf("-Floor needs an argument\n");
    return -1;
  }
  if (get_help || !filea || !fileb) {
    printf("\n");
    printf("usage: %s [options] file1.vtu file2.vtu\n", argv[0]);
    printf("    or:  %s [-help]             (usage)\n", argv[0]);
    printf("\n");
    printf("    -help (Print this summary and exit.)\n");
    printf("    -tolerance <$val> (Overrides the default tolerance of 1.0E-6.)\n");
    printf("    -Floor <$val> (Overrides the default floor tolerance of 0.0.)\n");
    osh_fini();
    return 0;
  }
  struct mesh* a = read_mesh_vtk(filea);
  struct mesh* b = read_mesh_vtk(fileb);
  unsigned differ = mesh_diff(a, b, tol, floor, allow_superset);
  free_mesh(a);
  free_mesh(b);
  osh_fini();
  if (differ)
    return 2;
  return 0;
}
