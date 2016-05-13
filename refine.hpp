#ifndef REFINE_HPP
#define REFINE_HPP

namespace omega_h {

struct mesh;

unsigned refine_by_size(struct mesh* m, double qual_floor);
void uniformly_refine(struct mesh* m, double qual_foor);

}

#endif
