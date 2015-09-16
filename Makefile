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

#these are source containing "library" functions,
#basically any source without a main() function
common_sources := \
star.c \
tables.c \
up_from_down.c \
ints.c \
vtk.c \
refine_topology.c \
splits_to_elements.c \
quality.c \
size.c \
bridge_graph.c \
refine_common.c \
refine_by_size.c \
concat.c \
indset.c \
measure_edges.c \
reflect_down.c \
refine_nodal.c \
refine_qualities.c \
doubles.c \
classify_box.c \
coarsen_by_size.c \
check_collapse_class.c \
coarsen_qualities.c \
coarsen_topology.c \
collapses_to_verts.c \
collapses_to_elements.c \
refine_class.c \
mesh.c \
graph.c \
field.c \
label.c \
warp_to_limit.c \
eval_field.c \
element_gradients.c \
jacobian.c \
recover_by_volume.c \
size_from_hessian.c \
subset.c \
adapt.c \
coarsen_common.c \
mark.c \
coarsen_slivers.c \
swap_slivers.c \
swap_common.c \
swap_qualities.c \
swap_topology.c \
edge_ring.c \
edge_swap.c \
loop.c

common_objects := $(common_sources:.c=.o)

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

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

#include the auto-generated dependency file for
#each source file
include $(depfiles)

.PHONY: default clean
