#ifndef MESH_H
#define MESH_H

#include "tag.h"

struct mesh;

struct const_up {
  unsigned const* const offsets;
  unsigned const* const adj;
  unsigned const* const directions;
};

struct mesh* new_mesh(unsigned elem_dim);
struct mesh* new_box_mesh(unsigned elem_dim);
void free_mesh(struct mesh* m);

unsigned mesh_dim(struct mesh* m);
unsigned mesh_count(struct mesh* m, unsigned dim);

void mesh_set_ents(struct mesh* m, unsigned dim, unsigned n, unsigned* verts);

unsigned const* mesh_ask_down(struct mesh* m,
    unsigned high_dim, unsigned low_dim);
struct const_up* mesh_ask_up(struct mesh* m,
    unsigned low_dim, unsigned high_dim);

struct const_graph* mesh_ask_star(struct mesh* m, unsigned low_dim,
    unsigned high_dim);
unsigned const* mesh_ask_dual(struct mesh* m);

struct const_tag* mesh_add_tag(struct mesh* m, unsigned dim, enum tag_type type,
    char const* name, unsigned ncomps, void* data);
void mesh_free_tag(struct mesh* m, unsigned dim, char const* name);
struct const_tag* mesh_find_tag(struct mesh* m, unsigned dim, char const* name);
unsigned mesh_count_tags(struct mesh* m, unsigned dim);
struct const_tag* mesh_get_tag(struct mesh* m, unsigned dim, unsigned i);

struct tags* mesh_tags(struct mesh* m, unsigned dim);
unsigned mesh_has_dim(struct mesh* m, unsigned dim);

struct parallel_mesh;

struct parallel_mesh* mesh_parallel(struct mesh* m);

#endif
