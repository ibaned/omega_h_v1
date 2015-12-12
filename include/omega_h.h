#ifndef OMEGA_H_H
#define OMEGA_H_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct osh_mesh* osh_t;

void osh_free(osh_t m);

osh_t osh_read_vtk(char const* filename);
void osh_write_vtk(osh_t m, char const* filename);

unsigned osh_dim(osh_t m);
unsigned osh_nelems(osh_t m);
unsigned osh_nverts(osh_t m);

unsigned const* osh_down(osh_t m, unsigned high_dim, unsigned low_dim);
unsigned const* osh_up(osh_t m, unsigned low_dim, unsigned high_dim);
unsigned const* osh_up_offs(osh_t m, unsigned low_dim, unsigned high_dim);
unsigned const* osh_up_dirs(osh_t m, unsigned low_dim, unsigned high_dim);

double const* osh_coords(osh_t m);

void osh_set_field(osh_t m, char const* name, unsigned ncomps, double* data);
void osh_new_field(osh_t m, char const* name, unsigned ncomps);
double* osh_get_field(osh_t m, char const* name);

void osh_mark_verts(osh_t m, unsigned class_dim, unsigned class_id,
    unsigned* marked);

void osh_add_label(osh_t m, char const* name, unsigned* data);
void osh_free_label(osh_t m, char const* name);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
