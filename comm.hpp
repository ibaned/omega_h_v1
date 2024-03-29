#ifndef COMM_HPP
#define COMM_HPP

namespace omega_h {

struct comm;

void comm_init(int* argc, char*** argv);
void comm_fini(void);

struct comm* comm_world(void);
struct comm* comm_self(void);

struct comm* comm_using(void);
void comm_use(struct comm* c);

struct comm* comm_split(struct comm* c, unsigned group, unsigned rank);

struct comm* comm_graph(struct comm* c,
    unsigned nout, unsigned const* out, unsigned const* outweights);
struct comm* comm_graph_exact(struct comm* c,
    unsigned nin, unsigned const* in, unsigned const* inweights,
    unsigned nout, unsigned const* out, unsigned const* outweights);
void comm_recvs(struct comm* c,
    unsigned* nin, unsigned** in, unsigned** inweights);
template <typename T>
void comm_exchange(struct comm* c,
    unsigned width,
    T const* out, unsigned const* outcounts, unsigned const* outoffsets,
    T* in, unsigned const* incounts, unsigned const* inoffsets);

void comm_sync_uint(struct comm* c, unsigned out, unsigned* in);
unsigned comm_bcast_uint(unsigned x);
void comm_bcast_chars(char* s, unsigned n);

void comm_free(struct comm* c);

unsigned comm_rank(void);
unsigned comm_size(void);

void comm_add_doubles(double* p, unsigned n);
double comm_add_double(double x);
double comm_max_double(double x);
double comm_min_double(double x);
unsigned long comm_add_ulong(unsigned long x);
unsigned long comm_exscan_ulong(unsigned long x);
unsigned long comm_max_ulong(unsigned long x);
unsigned comm_max_uint(unsigned x);

}

#endif
