#ifndef RV_MESH
#define RV_MESH

struct rv_mesh {
  unsigned elem_dim;
  unsigned nelems;
  unsigned nverts;
  unsigned* verts_of_elems;
  double* coords;
};

struct rv_mesh new_box_rv_mesh(unsigned dim);
void free_rv_mesh(struct rv_mesh m);

#endif
