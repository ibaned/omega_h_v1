CFLAGS = -Werror -Weverything -Wno-float-equal
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
bridge_graph.o
all: $(objects)

%: %.o $(common_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f *.dep* *.o

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
