#ifndef NODE_ELE_HPP
#define NODE_ELE_HPP

namespace omega_h {

struct mesh;

struct mesh* read_dot_node(char const* filename);
void read_dot_ele(struct mesh* m, char const* filename);
void write_dot_node(struct mesh* m, char const* filename);
void write_dot_ele(struct mesh* m, char const* filename);

}

#endif
