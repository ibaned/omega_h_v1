#config.mk is meant to contain the
#variables that commonly change,
#like CC and CFLAGS
#users edit config.mk but should never need
#to edit the Makefile
#this also separates the build configuration
#in the eye of the version control system,
#which can be told to ignore config.mk
include config.mk

#these are source containing "library" functions,
#basically any source without a main() function
sources := \
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
loop.c \
find_roots.c \
inertia.c \
derive_faces.c \
node_ele_io.c

ifeq "$(USE_MPI)" "yes"
sources += comm_mpi.c
else
sources += comm_serial.c
endif

targets ?= \
test_warp.exe \
test_warp_3d.exe \
vtk_surfer.exe

target_sources := $(patsubst %.exe,%.c,$(targets))
target_depfiles := $(patsubst %.exe,deps/%.dep,$(targets))

objects := $(patsubst %.c,objs/%.o,$(sources))
depfiles := $(patsubst %.c,deps/%.dep,$(sources))

#by default, the compilation target is to compile
#all .c files into objects
all: $(targets)

#general rule for an executable: link its object
#file with all the $(common_objects)
# $@ is the thing being built and $^ is all
#the things it depends on (the objects)
%.exe: objs/%.o $(objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

#cleanup removes dependency files, object files,
#and executables
clean:
	rm -rf deps/*.dep objs/*.o *.exe

#copied this mess from the GNU make documentation
#it generates dependency files from source files,
#and uses SED to change the rules
#such that the output is both an object file and a
#dependency file.
#it warrants further explanation:
#  cc -MM foo.c
#will produce a dependency line such as:
#  foo.o : foo.c foo.h bar.h
#the SED script changes this to:
#  objs/foo.o deps/foo.dep : foo.c foo.h bar.h
#$* is the same as the % in the rule pattern,
#"foo" in the example above.
#the $@ will insert the "deps/foo.dep"
deps/%.dep: %.c
	set -e; rm -f $@; \
	$(CPP) -MM $(CPPFLAGS) $< > $@.in; \
	sed 's,$*\.o,objs/$*.o $@,g' < $@.in > $@; \
	rm -f $@.in

#our rule for compiling a source file to an
#object, specifies that the object goes in objs/
objs/%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

#include the auto-generated dependency files for
#all source files.
#this is the funny recursion that keeps
#header file dependencies worked out at all times.
include $(depfiles) $(target_depfiles)

#"all" and "clean" are targets but not files or directories
.PHONY: all clean
