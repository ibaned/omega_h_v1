#ifndef MESH_H
#define MESH_H

struct graph;
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

struct const_field* mesh_add_field(struct mesh* m, unsigned dim, char const* name,
    unsigned ncomps, double* data);
void mesh_free_field(struct mesh* m, unsigned dim, char const* name);
struct const_field* mesh_find_field(struct mesh* m, unsigned dim,
    char const* name);
unsigned mesh_count_fields(struct mesh* m, unsigned dim);
struct const_field* mesh_get_field(struct mesh* m, unsigned dim, unsigned i);

struct const_label* mesh_add_label(struct mesh* m, unsigned dim, char const* name,
    unsigned* data);
struct const_label* mesh_find_label(struct mesh* m, unsigned dim,
    char const* name);
unsigned mesh_count_labels(struct mesh* m, unsigned dim);
struct const_label* mesh_get_label(struct mesh* m, unsigned dim, unsigned i);

#endif
