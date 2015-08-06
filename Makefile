CFLAGS = -O2 -std=c99 -Werror -Wall
LDLIBS = -lm

sources := $(wildcard *.c)
headers := $(wildcard *.h)
objects := $(patsubst %.c,%.o,$(sources))
depfiles := $(patsubst %.c,%.dep,$(sources))

common_objects := \
star.o \
tables.o \
up_from_down.o \
ints.o \
vtk.o \
refine_topology.o \
splits_to_elements.o \
quality.o \
size.o \
derive_edges.o \
bridge_graph.o \
refine_reduced.o \
concat.o \
indset.o \
measure_edges.o \
subset.o \
reflect_down.o \
refine_nodal.o \
rv_mesh.o \
refine_qualities.o

all: $(objects)

%.exe: %.o $(common_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f *.dep* *.o *.exe

#copied this mess from the GNU make documentation
#it generates dependency files from source files,
#and uses SED atrocities to change the rules
#such that output is both an object file and a
#dependency file
%.dep: %.c
	set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

#include the auto-generated dependency file for
#each source file
include $(depfiles)

.PHONY: default clean
