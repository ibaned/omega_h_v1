/* AUTO-GENERATED by print_swap_edges.exe !
   don't modify manually !

   tables describing, for each possible triangulation
   of an N-sided polygon, the interior edges of that
   triangulation.

   these will be used for full topology edge swapping
   to create intermediate entities in the interior of the
   cavity.

   we don't bother identifying unique edges because no caching
   is necessary for these intermediate entities */

#include "edge_swap.h"

LOOP_CONST static unsigned const edges_4_0[1 * 2] =
{2,0
};
LOOP_CONST static unsigned const edges_4_1[1 * 2] =
{1,3
};
LOOP_CONST static unsigned const edges_5_0[2 * 2] =
{2,0
,3,0
};
LOOP_CONST static unsigned const edges_5_1[2 * 2] =
{1,4
,1,3
};
LOOP_CONST static unsigned const edges_5_2[2 * 2] =
{2,0
,4,2
};
LOOP_CONST static unsigned const edges_5_3[2 * 2] =
{0,3
,3,1
};
LOOP_CONST static unsigned const edges_5_4[2 * 2] =
{1,4
,4,2
};
LOOP_CONST static unsigned const edges_6_0[3 * 2] =
{2,0
,3,0
,4,0
};
LOOP_CONST static unsigned const edges_6_1[3 * 2] =
{2,0
,2,5
,2,4
};
LOOP_CONST static unsigned const edges_6_2[3 * 2] =
{2,0
,3,0
,3,5
};
LOOP_CONST static unsigned const edges_6_3[3 * 2] =
{2,0
,0,4
,4,2
};
LOOP_CONST static unsigned const edges_6_4[3 * 2] =
{2,0
,2,5
,5,3
};
LOOP_CONST static unsigned const edges_6_5[3 * 2] =
{0,3
,4,0
,3,1
};
LOOP_CONST static unsigned const edges_6_6[3 * 2] =
{3,1
,1,5
,1,4
};
LOOP_CONST static unsigned const edges_6_7[3 * 2] =
{0,3
,3,5
,3,1
};
LOOP_CONST static unsigned const edges_6_8[3 * 2] =
{0,4
,3,1
,4,1
};
LOOP_CONST static unsigned const edges_6_9[3 * 2] =
{5,3
,3,1
,1,5
};
LOOP_CONST static unsigned const edges_6_10[3 * 2] =
{4,2
,1,5
,1,4
};
LOOP_CONST static unsigned const edges_6_11[3 * 2] =
{0,4
,4,2
,1,4
};
LOOP_CONST static unsigned const edges_6_12[3 * 2] =
{2,4
,5,2
,1,5
};
LOOP_CONST static unsigned const edges_6_13[3 * 2] =
{5,3
,5,2
,1,5
};
LOOP_CONST static unsigned const edges_7_0[4 * 2] =
{2,0
,3,0
,4,0
,5,0
};
LOOP_CONST static unsigned const edges_7_1[4 * 2] =
{2,0
,3,0
,3,6
,3,5
};
LOOP_CONST static unsigned const edges_7_2[4 * 2] =
{2,0
,3,0
,4,0
,4,6
};
LOOP_CONST static unsigned const edges_7_3[4 * 2] =
{2,0
,3,0
,0,5
,5,3
};
LOOP_CONST static unsigned const edges_7_4[4 * 2] =
{2,0
,3,0
,3,6
,6,4
};
LOOP_CONST static unsigned const edges_7_5[4 * 2] =
{2,0
,0,4
,5,0
,2,4
};
LOOP_CONST static unsigned const edges_7_6[4 * 2] =
{2,0
,4,2
,2,6
,2,5
};
LOOP_CONST static unsigned const edges_7_7[4 * 2] =
{2,0
,0,4
,4,6
,2,4
};
LOOP_CONST static unsigned const edges_7_8[4 * 2] =
{2,0
,0,5
,4,2
,5,2
};
LOOP_CONST static unsigned const edges_7_9[4 * 2] =
{2,0
,6,4
,4,2
,2,6
};
LOOP_CONST static unsigned const edges_7_10[4 * 2] =
{2,0
,5,3
,2,6
,2,5
};
LOOP_CONST static unsigned const edges_7_11[4 * 2] =
{2,0
,0,5
,5,3
,2,5
};
LOOP_CONST static unsigned const edges_7_12[4 * 2] =
{2,0
,3,5
,6,3
,2,6
};
LOOP_CONST static unsigned const edges_7_13[4 * 2] =
{2,0
,6,4
,6,3
,2,6
};
LOOP_CONST static unsigned const edges_7_14[4 * 2] =
{0,3
,4,0
,5,0
,1,3
};
LOOP_CONST static unsigned const edges_7_15[4 * 2] =
{0,3
,3,6
,3,5
,1,3
};
LOOP_CONST static unsigned const edges_7_16[4 * 2] =
{0,3
,4,0
,4,6
,1,3
};
LOOP_CONST static unsigned const edges_7_17[4 * 2] =
{0,5
,5,3
,0,3
,1,3
};
LOOP_CONST static unsigned const edges_7_18[4 * 2] =
{0,3
,3,6
,6,4
,1,3
};
LOOP_CONST static unsigned const edges_7_19[4 * 2] =
{0,4
,5,0
,3,1
,1,4
};
LOOP_CONST static unsigned const edges_7_20[4 * 2] =
{3,1
,4,1
,1,6
,1,5
};
LOOP_CONST static unsigned const edges_7_21[4 * 2] =
{0,4
,4,6
,3,1
,1,4
};
LOOP_CONST static unsigned const edges_7_22[4 * 2] =
{0,5
,3,1
,4,1
,5,1
};
LOOP_CONST static unsigned const edges_7_23[4 * 2] =
{6,4
,3,1
,4,1
,1,6
};
LOOP_CONST static unsigned const edges_7_24[4 * 2] =
{5,3
,3,1
,1,6
,1,5
};
LOOP_CONST static unsigned const edges_7_25[4 * 2] =
{0,5
,5,3
,3,1
,1,5
};
LOOP_CONST static unsigned const edges_7_26[4 * 2] =
{3,5
,6,3
,3,1
,1,6
};
LOOP_CONST static unsigned const edges_7_27[4 * 2] =
{6,4
,6,3
,3,1
,1,6
};
LOOP_CONST static unsigned const edges_7_28[4 * 2] =
{0,4
,5,0
,4,2
,1,4
};
LOOP_CONST static unsigned const edges_7_29[4 * 2] =
{4,2
,1,6
,1,5
,1,4
};
LOOP_CONST static unsigned const edges_7_30[4 * 2] =
{0,4
,4,6
,4,2
,1,4
};
LOOP_CONST static unsigned const edges_7_31[4 * 2] =
{0,5
,4,2
,1,4
,5,1
};
LOOP_CONST static unsigned const edges_7_32[4 * 2] =
{6,4
,4,2
,1,6
,1,4
};
LOOP_CONST static unsigned const edges_7_33[4 * 2] =
{4,2
,5,2
,1,6
,1,5
};
LOOP_CONST static unsigned const edges_7_34[4 * 2] =
{0,5
,4,2
,5,2
,1,5
};
LOOP_CONST static unsigned const edges_7_35[4 * 2] =
{4,2
,2,5
,6,2
,1,6
};
LOOP_CONST static unsigned const edges_7_36[4 * 2] =
{6,4
,4,2
,6,2
,1,6
};
LOOP_CONST static unsigned const edges_7_37[4 * 2] =
{5,3
,5,2
,1,6
,1,5
};
LOOP_CONST static unsigned const edges_7_38[4 * 2] =
{0,5
,5,3
,5,2
,1,5
};
LOOP_CONST static unsigned const edges_7_39[4 * 2] =
{5,3
,2,5
,6,2
,1,6
};
LOOP_CONST static unsigned const edges_7_40[4 * 2] =
{3,5
,6,3
,6,2
,1,6
};
LOOP_CONST static unsigned const edges_7_41[4 * 2] =
{6,4
,6,3
,6,2
,1,6
};
LOOP_CONST static unsigned const* const edges_4[2] =
{edges_4_0
,edges_4_1
};
LOOP_CONST static unsigned const* const edges_5[5] =
{edges_5_0
,edges_5_1
,edges_5_2
,edges_5_3
,edges_5_4
};
LOOP_CONST static unsigned const* const edges_6[14] =
{edges_6_0
,edges_6_1
,edges_6_2
,edges_6_3
,edges_6_4
,edges_6_5
,edges_6_6
,edges_6_7
,edges_6_8
,edges_6_9
,edges_6_10
,edges_6_11
,edges_6_12
,edges_6_13
};
LOOP_CONST static unsigned const* const edges_7[42] =
{edges_7_0
,edges_7_1
,edges_7_2
,edges_7_3
,edges_7_4
,edges_7_5
,edges_7_6
,edges_7_7
,edges_7_8
,edges_7_9
,edges_7_10
,edges_7_11
,edges_7_12
,edges_7_13
,edges_7_14
,edges_7_15
,edges_7_16
,edges_7_17
,edges_7_18
,edges_7_19
,edges_7_20
,edges_7_21
,edges_7_22
,edges_7_23
,edges_7_24
,edges_7_25
,edges_7_26
,edges_7_27
,edges_7_28
,edges_7_29
,edges_7_30
,edges_7_31
,edges_7_32
,edges_7_33
,edges_7_34
,edges_7_35
,edges_7_36
,edges_7_37
,edges_7_38
,edges_7_39
,edges_7_40
,edges_7_41
};
LOOP_CONST unsigned const* const* const swap_int_edges[MAX_EDGE_SWAP + 1] =
{0
,0
,0
,0
,edges_4
,edges_5
,edges_6
,edges_7
};
LOOP_CONST unsigned const swap_nint_edges[MAX_EDGE_SWAP + 1] =
{0
,0
,0
,0
,1
,2
,3
,4
};
