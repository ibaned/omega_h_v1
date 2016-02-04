#include "tables.h"

#include "arrays.h"
#include "loop.h"

static unsigned const box_1d_conn[1 * 2] = {
  0, 1
};

static double const box_1d_coord[2 * 3] = {
  0, 0, 0,
  1, 0, 0
};

static unsigned const box_2d_conn[2 * 3] = {
  0, 1, 2,
  2, 3, 0
};

static double const box_2d_coord[4 * 3] = {
  0, 0, 0,
  1, 0, 0,
  1, 1, 0,
  0, 1, 0
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

unsigned const the_box_nelems[4] = {
  1,
  1,
  2,
  6
};

unsigned const the_box_nverts[4] = {
  1,
  2,
  4,
  8
};

unsigned const the_down_degrees[4][4] = {
  { 1, 0, 0, 0},
  { 2, 1, 0, 0},
  { 3, 3, 1, 0},
  { 4, 6, 4, 1},
};

/* beyond this comment lie the vast tables
 * in which the Programmer in His wisdom
 * encoded the local Canonical Orders for all
 * simplicial elements and their boundaries
 */
static unsigned const vvv0[1] = {0};
static unsigned const* const vvv[1] = {vvv0};
static unsigned const* const* const vv_[1] = {vvv};
static unsigned const* const* const* const v__[1] = {vv_};
static unsigned const evv0[1] = {0};
static unsigned const evv1[2] = {1};
static unsigned const* const evv[2] = {evv0,evv1};
static unsigned const* const* const ev_[1] = {evv};
static unsigned const eev0[2] = {0,1};
static unsigned const* const eev[1] = {eev0};
static unsigned const eee0[1] = {0};
static unsigned const* const eee[1] = {eee0};
static unsigned const* const* const ee_[2] = {eev,eee};
static unsigned const* const* const* const e__[2] = {ev_,ee_};
static unsigned const fvv0[1] = {0};
static unsigned const fvv1[1] = {1};
static unsigned const fvv2[1] = {2};
static unsigned const* const fvv[3] = {fvv0,fvv1,fvv2};
static unsigned const* const* const fv_[1] = {fvv};
static unsigned const fev0[2] = {0,1}; /* cyclic */
static unsigned const fev1[2] = {1,2};
static unsigned const fev2[2] = {2,0};
static unsigned const* const fev[3] = {fev0,fev1,fev2};
static unsigned const fee0[1] = {0};
static unsigned const fee1[1] = {1};
static unsigned const fee2[1] = {2};
static unsigned const* const fee[3] = {fee0,fee1,fee2};
static unsigned const* const* const fe_[2] = {fev,fee};
static unsigned const ffv0[3] = {0,1,2};
static unsigned const* const ffv[1] = {ffv0};
static unsigned const ffe0[3] = {0,1,2};
static unsigned const* const ffe[1] = {ffe0};
static unsigned const fff0[1] = {0};
static unsigned const* const fff[1] = {fff0};
static unsigned const* const* const ff_[3] = {ffv,ffe,fff};
static unsigned const* const* const* const f__[3] = {fv_,fe_,ff_};
static unsigned const rvv0[1] = {0};
static unsigned const rvv1[1] = {1};
static unsigned const rvv2[1] = {2};
static unsigned const rvv3[1] = {3};
static unsigned const* const rvv[4] = {rvv0,rvv1,rvv2,rvv3};
static unsigned const* const* const rv_[1] = {rvv};
static unsigned const rev0[2] = {0,1};
static unsigned const rev1[2] = {1,2};
static unsigned const rev2[2] = {2,0};
static unsigned const rev3[2] = {0,3};
static unsigned const rev4[2] = {1,3};
static unsigned const rev5[2] = {2,3};
static unsigned const* const rev[6] = {rev0,rev1,rev2,rev3,rev4,rev5};
static unsigned const ree0[1] = {0};
static unsigned const ree1[1] = {1};
static unsigned const ree2[1] = {2};
static unsigned const ree3[1] = {3};
static unsigned const ree4[1] = {4};
static unsigned const ree5[1] = {5};
static unsigned const* const ree[6] = {ree0,ree1,ree2,ree3,ree4,ree5};
static unsigned const* const* const re_[2] = {rev,ree};
static unsigned const rfv0[3] = {0,1,2};
static unsigned const rfv1[3] = {1,0,3}; /* all curl inward */
static unsigned const rfv2[3] = {2,1,3};
static unsigned const rfv3[3] = {0,2,3};
static unsigned const* const rfv[4] = {rfv0,rfv1,rfv2,rfv3};
static unsigned const rfe0[3] = {0,1,2};
static unsigned const rfe1[3] = {0,3,4}; /* all curl inward */
static unsigned const rfe2[3] = {1,4,5};
static unsigned const rfe3[3] = {2,5,3};
static unsigned const* const rfe[4] = {rfe0,rfe1,rfe2,rfe3};
static unsigned const rff0[1] = {0};
static unsigned const rff1[1] = {1};
static unsigned const rff2[1] = {2};
static unsigned const rff3[1] = {3};
static unsigned const* const rff[4] = {rff0,rff1,rff2,rff3};
static unsigned const* const* const rf_[3] = {rfv,rfe,rff};
static unsigned const rrv0[4] = {0,1,2,3};
static unsigned const* const rrv[1] = {rrv0};
static unsigned const rre0[6] = {0,1,2,3,4,5};
static unsigned const* const rre[1] = {rre0};
static unsigned const rrf0[4] = {0,1,2,3};
static unsigned const* const rrf[1] = {rrf0};
static unsigned const rrr0[1] = {0};
static unsigned const* const rrr[1] = {rrr0};
static unsigned const* const* const rr_[4] = {rrv,rre,rrf,rrr};
static unsigned const* const* const* const r__[4] = {rv_,re_,rf_,rr_};
unsigned const* const* const* const* const the_canonical_orders[4] = {
  v__,
  e__,
  f__,
  r__
};

static unsigned const oev[2] = {1,0};
static unsigned const* const oe_[1] = {oev};
static unsigned const ofv[3] = {1,2,0};
static unsigned const ofe[3] = {2,0,1};
static unsigned const* const of_[2] = {ofv,ofe};
static unsigned const orv[4] = {2,3,1,0};
static unsigned const ore[6] = {5,3,4,1,2,0};
static unsigned const orf[4] = {3,2,0,1};
static unsigned const* const or_[3] = {orv,ore,orf};
unsigned const* const* const the_opposite_orders[4] = {
  0,
  oe_,
  of_,
  or_
};

unsigned** orders_to_device(unsigned dim1, unsigned dim2, unsigned dim3)
{
  unsigned deg1 = the_down_degrees[dim1][dim2];
  unsigned deg2 = the_down_degrees[dim2][dim3];
  unsigned* a[MAX_DOWN];
  for (unsigned i = 0; i < deg1; ++i)
    a[i] = uints_to_device(
        the_canonical_orders[dim1][dim2][dim3][i], deg2);
  return uintptrs_to_device(a, deg1);
}

void free_orders(unsigned** a, unsigned dim1, unsigned dim2)
{
  unsigned deg1 = the_down_degrees[dim1][dim2];
  unsigned** b = LOOP_TO_HOST(unsigned*, a, deg1);
  loop_free(a);
  for (unsigned i = 0; i < deg1; ++i)
    loop_free(b[i]);
  loop_host_free(b);
}

static char const* const the_ent_names[4] = {
  "vertex",
  "edge",
  "triangle",
  "tet"
};

static char const* const the_plural_ent_names[4] = {
  "vertices",
  "edges",
  "triangles",
  "tets"
};

char const* get_ent_name(unsigned dim, unsigned long n)
{
  if (n == 1)
    return the_ent_names[dim];
  return the_plural_ent_names[dim];
}
