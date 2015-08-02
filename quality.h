#ifndef QUALITY_H
#define QUALITY_H

typedef double (*quality_function)(double coords[][3]);

extern quality_function const the_quality_functions[4];

#endif
