#ifndef GLOBAL_HPP
#define GLOBAL_HPP

namespace omega_h {

struct mesh;

unsigned long* globalize_offsets(unsigned* local, unsigned n);
void global_to_linpart(unsigned long const* global, unsigned n,
    unsigned long total, unsigned nparts,
    unsigned** p_part, unsigned** p_local);
unsigned linpart_size(unsigned long total, unsigned nparts, unsigned part);

}

#endif
