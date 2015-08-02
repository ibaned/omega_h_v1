CFLAGS = -Weverything
sources := $(wildcard *.c)
headers := $(wildcard *.h)
objects := $(patsubst %.c,%.o,$(sources))
depfiles := $(patsubst %.c,%.d,$(sources))
common_objects := \
star.o \
tables.o \
up_from_down.o \
ints.o \
vtk.o \
refine_topology.o \
splits_to_elements.o
all: $(objects)
test_up_from_down: test_up_from_down.o $(common_objects)
test_vtk: test_vtk.o $(common_objects)
test_refine_topology: test_refine_topology.o $(common_objects)
clean:
	rm -f $(objects) $(depfiles)

#copied this mess from the GNU make documentation
#it generates dependency files from source files,
#and uses SED atrocities to change the rules
#such that output is both an object file and a
#dependency file
%.d: %.c
	set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

#include the auto-generated dependency file for
#each source file
include $(depfiles)

.PHONY: default clean
