#ifndef COLLAPSE_H
#define COLLAPSE_H

enum {
  DONT_COLLAPSE  = 0,
  COLLAPSE_LEFT  = 1,
  COLLAPSE_RIGHT = 2,
  COLLAPSE_BOTH  = 3,
};

static inline unsigned collapses(unsigned code, unsigned dir)
{
  return (code & (((unsigned)1) << dir)) != 0;
}

static inline unsigned dont_collapse(unsigned code, unsigned dir)
{
  return code & ~(((unsigned)1) << dir);
}

#endif
