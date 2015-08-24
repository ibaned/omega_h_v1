CFLAGS = -g -O2 -std=c99 -Werror -Wall
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
bridge_graph.o \
refine_common.o \
refine_by_size.o \
concat.o \
indset.o \
measure_edges.o \
reflect_down.o \
refine_nodal.o \
refine_qualities.o \
doubles.o \
classify_box.o \
coarsen_by_size.o \
check_collapse_class.o \
coarsen_qualities.o \
coarsen_topology.o \
collapses_to_verts.o \
collapses_to_elements.o \
verify.o \
split_sliver_tris.o \
bad_elem_keys.o \
collect_keys.o \
refine_class.o \
derive_faces.o \
mesh.o \
graph.o \
field.o \
label.o \
inside.o \
warp_to_limit.o

all: $(objects)

%.exe: %.o $(common_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_refine_by_size.exe: test_refine_by_size.o $(common_objects)
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
