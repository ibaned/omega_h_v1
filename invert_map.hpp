#ifndef INVERT_MAP_HPP
#define INVERT_MAP_HPP

namespace omega_h {

/* given an array (in) which represents a
   map from the integers A=[0,nin) to the integers B=[0,nout),
   as defined by b = in[a], this function creates the inverse
   map as two arrays (out) and (offsets), where
   {out[offsets[b]], out[offsets[b]+1], ..., out[offsets[b+1]-1]}
   are the integers (a) such that in[a]=b
*/

void invert_map(
    unsigned nin,
    unsigned const* in,
    unsigned nout,
    unsigned** p_out,
    unsigned** p_offsets);

}

#endif
