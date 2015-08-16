#ifndef MESH_H
#define MESH_H

#include "graph.h"
#include "field.h"
#include "label.h"

struct up {
  unsigned* offsets;
  unsigned* adj;
  unsigned* directions;
};

struct up* new_up(unsigned* offsets, unsigned* adj, unsigned* directions);
void free_up(struct up* u);

struct mesh* new_mesh(unsigned elem_dim);
struct mesh* new_box_mesh(unsigned elem_dim);
unsigned mesh_dim(struct mesh* m);
unsigned mesh_count(struct mesh* m, unsigned dim);
void free_mesh(struct mesh* m);
struct field* mesh_find_nodal_field(struct mesh* m, char const* name);
unsigned* mesh_ask_down(struct mesh* m, unsigned high_dim, unsigned low_dim);
struct up* mesh_ask_up(struct mesh* m, unsigned low_dim, unsigned high_dim);
struct graph* mesh_ask_star(struct mesh* m, unsigned low_dim, unsigned high_dim);
unsigned* mesh_ask_dual(struct mesh* m);
void mesh_set_down(struct mesh* m, unsigned high_dim, unsigned low_dim,
    unsigned* adj);
void mesh_set_up(struct mesh* m, unsigned low_dim, unsigned high_dim,
    struct up* adj);
void mesh_set_star(struct mesh* m, unsigned low_dim, unsigned high_dim,
    struct graph* adj);
void mesh_set_dual(struct mesh* m, unsigned* adj);
void mesh_free_down(struct mesh* m, unsigned high_dim, unsigned low_dim);
void mesh_free_up(struct mesh* m, unsigned low_dim, unsigned high_dim);
void mesh_free_star(struct mesh* m, unsigned low_dim, unsigned high_dim);
void mesh_free_dual(struct mesh* m);
void mesh_set_ents(struct mesh* m, unsigned dim, unsigned n, unsigned* verts);
struct field* mesh_add_nodal_field(struct mesh* m, char const* name,
    unsigned ncomps);

#endif
