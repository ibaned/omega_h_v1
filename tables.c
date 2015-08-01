#include "tables.h"

static unsigned const box_1d_conn[1 * 2] = {
  0, 1
};

static double const box_1d_coord[2 * 1] = {
  0,
  1
};

static unsigned const box_2d_conn[2 * 3] = {
  0, 1, 2,
  2, 3, 0
};

static double const box_2d_coord[4 * 2] = {
  0, 0,
  1, 0,
  1, 1,
  0, 1
};

static unsigned const box_3d_conn[6 * 4] = {
  0, 1, 2, 6,
  2, 3, 0, 6,
  0, 3, 7, 6,
  7, 4, 0, 6,
  0, 4, 5, 6,
  5, 1, 0, 6
};

static double const box_3d_coord[8 * 3] = {
  0, 0, 0,
  1, 0, 0,
  1, 1, 0,
  0, 1, 0,
  0, 0, 1,
  1, 0, 1,
  1, 1, 1,
  0, 1, 1
};

unsigned const* const the_box_conns[4] = {
  0,
  box_1d_conn,
  box_2d_conn,
  box_3d_conn
};

double const* const the_box_coords[4] = {
  0,
  box_1d_coord,
  box_2d_coord,
  box_3d_coord
};

unsigned const the_down_degrees[4][4] = {
  { 1, 0, 0, 0},
  { 2, 1, 0, 0},
  { 3, 3, 1, 0},
  { 4, 6, 4, 1},
};
