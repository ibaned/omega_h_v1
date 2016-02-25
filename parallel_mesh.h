#ifndef PARALLEL_MESH_H
#define PARALLEL_MESH_H

struct mesh;
struct parallel_mesh;
struct exchanger;

struct parallel_mesh* new_parallel_mesh(void);
void free_parallel_mesh(struct parallel_mesh* m);

unsigned long const* mesh_ask_globals(struct mesh* m, unsigned dim);
unsigned const* mesh_ask_own_ranks(struct mesh* m, unsigned dim);
unsigned const* mesh_ask_own_ids(struct mesh* m, unsigned dim);
struct exchanger* mesh_ask_exchanger(struct mesh* m, unsigned dim);

void mesh_global_renumber(struct mesh* m, unsigned dim);

void mesh_conform_doubles(struct mesh* m, unsigned dim, unsigned width,
    double** a);
void mesh_conform_uints(struct mesh* m, unsigned dim, unsigned width,
    unsigned** a);
void mesh_conform_ulongs(struct mesh* m, unsigned dim, unsigned width,
    unsigned long** a);

void mesh_conform_tag(struct mesh* m, unsigned dim, const char* name);
void mesh_accumulate_tag(struct mesh* m, unsigned dim, const char* name);

unsigned mesh_ghost_layers(struct mesh* m);

void mesh_set_ghost_layers(struct mesh* m, unsigned n);

void mesh_set_globals(struct mesh* m, unsigned dim, unsigned long* new_globals);
void mesh_set_own_ranks(struct mesh* m, unsigned dim, unsigned* new_owners);

void mesh_tag_globals(struct mesh* m, unsigned dim);
void mesh_tag_own_rank(struct mesh* m, unsigned dim);
void mesh_tag_own_id(struct mesh* m, unsigned dim);

void mesh_parallel_to_tags(struct mesh* m, unsigned dim);
void mesh_parallel_untag(struct mesh* m, unsigned dim);
void mesh_parallel_from_tags(struct mesh* m, unsigned dim);

void mesh_partition_out(struct mesh** p_m, unsigned factor);
void mesh_partition_in(struct mesh** p_m, unsigned factor);

struct mesh* read_and_partition_serial_mesh(char const* filename);

unsigned* mesh_get_owned(struct mesh* m, unsigned dim);

void mesh_reduce_doubles_max(struct mesh* m, unsigned dim, unsigned width,
    double** a);

#endif
