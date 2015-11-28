#ifndef PARALLEL_MESH_H
#define PARALLEL_MESH_H

struct mesh;
struct parallel_mesh;
struct exchanger;

struct parallel_mesh* new_parallel_mesh(struct mesh* m);
void free_parallel_mesh(struct parallel_mesh* m);

unsigned long const* mesh_ask_global(struct mesh* m, unsigned dim);
unsigned const* mesh_ask_own_ranks(struct mesh* m, unsigned dim);
unsigned const* mesh_ask_own_ids(struct mesh* m, unsigned dim);
struct exchanger* mesh_ask_exchanger(struct mesh* m, unsigned dim);

unsigned long* conform_ulongs(struct exchanger* ex, unsigned width,
    unsigned long const* a);
double* conform_doubles(struct exchanger* ex, unsigned width,
    double const* a);
void mesh_global_renumber(struct mesh* m, unsigned dim);

#endif
