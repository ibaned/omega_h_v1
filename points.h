#ifndef POINTS_H
#define POINTS_H

struct points {
  unsigned* offsets;
  unsigned* adj;
  double* coords;
};

struct const_points {
  unsigned const* const offsets;
  unsigned const* const adj;
  double const* const coords;
};

struct points* new_points(unsigned* offsets, unsigned* adj, double* coords);
void free_points(struct points* ps);

#endif
