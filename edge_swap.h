#ifndef EDGE_SWAP_H
#define EDGE_SWAP_H

#define MAX_EDGE_SWAP 7

struct swap_choice {
  unsigned code;
  double quality;
};

struct swap_choice choose_edge_swap(
    unsigned ring_size,
    double edge_x[2][3],
    double ring_x[][3],
    double good_qual,
    double valid_qual);

#endif
