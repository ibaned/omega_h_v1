#ifndef OMEGA_H_H
#define OMEGA_H_H

typedef struct osh_mesh* osh_t;

osh_t osh_read_vtk(char const* filename);
unsigned const* osh_ask_down(osh_t m, unsigned high_dim, unsigned low_dim);
double const* osh_ask_coords(osh_t m);
void osh_free(osh_t m);

#endif
