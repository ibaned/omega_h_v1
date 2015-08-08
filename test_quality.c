#include "quality.h"
#include <math.h>
#include <stdio.h>

int main()
{
  double tri_coords[3][3] = {
    {-1,0,0},
    {0,sqrt(3),0},
    {1,0,0}
  };
  printf("regular triangle quality = %f\n", triangle_quality(tri_coords));
  double tet_coords[4][3] = {
    { 1,-1,-1},
    { 1, 1, 1},
    {-1, 1,-1},
    {-1,-1, 1} 
  };
  printf("regular tet quality = %f\n", tet_quality(tet_coords));
}
