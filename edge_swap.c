#include "edge_swap.h"
#include "algebra.h"
#include "quality.h"
#include "tables.h"
#include <assert.h>

unsigned const swap_mesh_sizes[MAX_EDGE_SWAP+1] =
{0 //0
,0 //1
,0 //2
,1 //3
,2 //4
,3 //5
,4 //6
,5 //7
};

static unsigned const mesh_counts[MAX_EDGE_SWAP+1] =
{0  //0
,0  //1
,0  //2
,1  //3
,2  //4
,5  //5
,14 //6
,42 //7
};

static unsigned const unique_triangle_counts[MAX_EDGE_SWAP+1] =
{0  //0
,0  //1
,0  //2
,1  //3
,4  //4
,10 //5
,20 //6
,35 //7
};

#define MAX_UNIQUE_TRIS 35

static unsigned const triangles_3[1][3] = {{0,1,2}};

static unsigned const meshes_3[1] =
{0
};

static unsigned const triangles_4[4][3] =
{{0,1,2}
,{0,2,3}
,{0,1,3}
,{1,2,3}
};

static unsigned const meshes_4[2 * 2] =
{0,1
,2,3
};

static unsigned const triangles_5[10][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,1,4}
,{1,3,4}
,{1,2,3}
,{2,3,4}
,{0,2,4}
,{0,1,3}
,{1,2,4}
};

static unsigned const meshes_5[5 * 3] =
{0,1,2
,3,4,5
,0,6,7
,2,5,8
,3,6,9
};

static unsigned const triangles_6[20][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,4,5}
,{0,2,5}
,{2,4,5}
,{2,3,4}
,{0,3,5}
,{3,4,5}
,{0,2,4}
,{2,3,5}
,{1,2,3}
,{0,1,3}
,{0,1,5}
,{1,4,5}
,{1,3,4}
,{0,1,4}
,{1,3,5}
,{1,2,4}
,{1,2,5}
};

static unsigned const meshes_6[14 * 4] =
{0, 1 ,2 ,3
,0, 4 ,5 ,6
,0, 1 ,7 ,8
,0, 3 ,6 ,9
,0, 4 ,8 ,10
,2, 3 ,11,12
,11,13,14,15
,7 ,8 ,11,12
,3 ,11,15,16
,8 ,11,13,17
,6 ,13,14,18
,3 ,6 ,16,18
,5 ,6 ,13,19
,8 ,10,13,19
};

static unsigned const triangles_7[35][3] =
{{0,1,2}
,{0,2,3}
,{0,3,4}
,{0,4,5}
,{0,5,6}
,{0,3,6}
,{3,5,6}
,{3,4,5}
,{0,4,6}
,{4,5,6}
,{0,3,5}
,{3,4,6}
,{0,2,4}
,{2,3,4}
,{0,2,6}
,{2,5,6}
,{2,4,5}
,{0,2,5}
,{2,4,6}
,{2,3,5}
,{2,3,6}
,{0,1,3}
,{1,2,3}
,{0,1,4}
,{1,3,4}
,{0,1,6}
,{1,5,6}
,{1,4,5}
,{0,1,5}
,{1,4,6}
,{1,3,5}
,{1,3,6}
,{1,2,4}
,{1,2,5}
,{1,2,6}
};

static unsigned const meshes_7[42 * 5] =
{0 ,1 ,2 ,3 ,4
,0 ,1 ,5 ,6 ,7
,0 ,1 ,2 ,8 ,9
,0 ,1 ,4 ,7 ,10
,0 ,1 ,5 ,9 ,11
,0 ,3 ,4 ,12,13
,0 ,13,14,15,16
,0 ,8 ,9 ,12,13
,0 ,4 ,13,16,17
,0 ,9 ,13,14,18
,0 ,7 ,14,15,19
,0 ,4 ,7 ,17,19
,0 ,6 ,7 ,14,20
,0 ,9 ,11,14,20
,2 ,3 ,4 ,21,22
,5 ,6 ,7 ,21,22
,2 ,8 ,9 ,21,22
,4 ,7 ,10,21,22
,5 ,9 ,11,21,22
,3 ,4 ,22,23,24
,22,24,25,26,27
,8 ,9 ,22,23,24
,4 ,22,24,27,28
,9 ,22,24,25,29
,7 ,22,25,26,30
,4 ,7 ,22,28,30
,6 ,7 ,22,25,31
,9 ,11,22,25,31
,3 ,4 ,13,23,32
,13,25,26,27,32
,8 ,9 ,13,23,32
,4 ,13,27,28,32
,9 ,13,25,29,32
,13,16,25,26,33
,4 ,13,16,28,33
,13,15,16,25,34
,9 ,13,18,25,34
,7 ,19,25,26,33
,4 ,7 ,19,28,33
,7 ,15,19,25,34
,6 ,7 ,20,25,34
,9 ,11,20,25,34
};

/* array [8] of pointer to array [3] of unsigned const */
typedef unsigned const tri_t[3];
static tri_t const* const triangles[MAX_EDGE_SWAP+1] =
{0
,0
,0
,triangles_3
,triangles_4
,triangles_5
,triangles_6
,triangles_7
};

static unsigned const* const meshes[MAX_EDGE_SWAP+1] =
{0
,0
,0
,meshes_3
,meshes_4
,meshes_5
,meshes_6
,meshes_7
};

struct swap_choice choose_edge_swap(
    unsigned ring_size,
    double edge_x[2][3],
    double ring_x[][3],
    double good_qual)
{
  unsigned tris_per_mesh = swap_mesh_sizes[ring_size];
  unsigned nmeshes = mesh_counts[ring_size];
  unsigned ntris = unique_triangle_counts[ring_size];
  unsigned const* mesh = meshes[ring_size];
  tri_t const* tris = triangles[ring_size];
  unsigned char cached[MAX_UNIQUE_TRIS];
  double cache[MAX_UNIQUE_TRIS];
  for (unsigned i = 0; i < ntris; ++i)
    cached[i] = 0;
  struct swap_choice out;
  out.code = INVALID;
  out.padding__ = 0;
  out.quality = -1;
  for (unsigned i = 0; i < nmeshes; ++i) {
    double mesh_minq = 1;
    for (unsigned j = 0; j < tris_per_mesh; ++j) {
      unsigned tri = mesh[j];
      double tri_minq;
      if (!cached[tri]) {
        unsigned const* tri_verts = tris[tri];
        double tet_x[4][3];
        for (unsigned k = 0; k < 3; ++k) {
          unsigned vert = tri_verts[k];
          copy_vector(ring_x[vert], tet_x[k], 3);
        }
        tri_minq = 1;
        for (unsigned k = 0; k < 2; ++k) {
          copy_vector(edge_x[1 - k], tet_x[3], 3);
          double tet_q = tet_quality(tet_x);
          if (tet_q < tri_minq)
            tri_minq = tet_q;
          swap_vectors(tet_x[0], tet_x[1], 3);
        }
        cached[tri] = 1;
        cache[tri] = tri_minq;
      }
      tri_minq = cache[tri];
      if (tri_minq < mesh_minq)
        mesh_minq = tri_minq;
      if (mesh_minq < 0)
        break;
    }
    if (mesh_minq < 0)
      /* FIXME: the following is a BUG */
      continue;
    if (mesh_minq > out.quality) {
      out.code = i;
      out.quality = mesh_minq;
    }
    /* FIXME: we may not want to do this */
    if (out.quality >= good_qual)
      break;
    mesh += tris_per_mesh;
  }
  return out;
}

void apply_edge_swap(
    unsigned ring_size,
    unsigned code,
    unsigned const edge_v[2],
    unsigned const ring_v[],
    unsigned out[])
{
  unsigned* p = out;
  unsigned tris_per_mesh = swap_mesh_sizes[ring_size];
  unsigned const* mesh = meshes[ring_size] + code * tris_per_mesh;
  tri_t const* tris = triangles[ring_size];
  for (unsigned i = 0; i < tris_per_mesh; ++i) {
    unsigned tet_verts[4];
    unsigned tri = mesh[i];
    unsigned const* tri_verts = tris[tri];
    for (unsigned j = 0; j < 3; ++j)
      tet_verts[j] = ring_v[tri_verts[j]];
    for (unsigned j = 0; j < 2; ++j) {
      tet_verts[3] = edge_v[1 - j];
      for (unsigned k = 0; k < 4; ++k)
        *out++ = tet_verts[k];
      unsigned tmp = tet_verts[0];
      tet_verts[0] = tet_verts[1];
      tet_verts[1] = tmp;
    }
  }
  unsigned ntets = tris_per_mesh * 2;
  assert(out - p == ntets * 4);
}
