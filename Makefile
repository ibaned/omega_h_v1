#config.mk is meant to contain the
#variables that commonly change,
#like CC and CFLAGS
#users edit config.mk but should never need
#to edit the Makefile
#this also separates the build configuration
#in the eye of the version control system,
#which can be told to ignore config.mk
include config.mk

#take *all* files ending in .c, use that as
#the list of sources to compile
sources := $(wildcard *.c)
#the list of objects to compile is derived
#from the list of source by changing .c to .o
objects := $(sources:.c=.o)
#the list of dependency files is also derived
#from the list of source by changing .c to .dep
depfiles := $(sources:.c=.dep)

#these are objects containing "library" functions,
#basically any object without a main() function
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
refine_class.o \
mesh.o \
graph.o \
field.o \
label.o \
warp_to_limit.o \
eval_field.o \
element_gradients.o \
jacobian.o \
recover_by_volume.o \
size_from_hessian.o \
subset.o \
adapt.o \
coarsen_common.o \
mark.o \
coarsen_slivers.o \
swap_slivers.o \
swap_common.o \
swap_qualities.o \
swap_topology.o \
edge_ring.o \
edge_swap.o \
loop.o

#by default, the compilation target is to compile
#all .c files into objects
all: $(objects)

#general rule for an executable: link its object
#file with all the $(common_objects)
# $@ is the thing being built and $^ is all
#the things it depends on (the objects)
%.exe: %.o $(common_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

#cleanup removes dependency files, object files,
#and executables
clean:
	rm -f *.dep* *.o *.exe

#copied this mess from the GNU make documentation
#it generates dependency files from source files,
#and uses SED atrocities to change the rules
#such that output is both an object file and a
#dependency file
#it warrants further explanation:
#  cc -MM foo.c
#will produce a dependency line such as:
#  foo.o : foo.c foo.h bar.h
#the SED script changes this to:
#  foo.o foo.dep : foo.c foo.h bar.h
#$* is the same as the % in the rule pattern,
#"foo" in the example above.
#\($*\)\.o matches foo.o and stores "foo" for later
#[ :]* matches all the space and the colon
#between foo.o and foo.c
#\1.o will print out the "foo" that was stored
#from before and $@ will insert the "foo.dep"
#$$$$ gets replaced with the PID for `make`
#to create a temporary file...
%.dep: %.c
	set -e; rm -f $@; \
	$(CPP) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

#include the auto-generated dependency file for
#each source file
include $(depfiles)

.PHONY: default clean
