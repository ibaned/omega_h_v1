#The user configuration is placed in config.mk
#Users edit config.mk but should never need
#to edit the Makefile.
#This also separates the build configuration
#in the eye of the version control system,
#which does not keep track of config.mk
include config.mk

test_sources := \
test_coarsen_3d.c \
test_coarsen_by_size.c \
test_eigenval.c \
test_form_cloud.c \
test_grad.c \
test_inertia.c \
test_node_ele.c \
test_quality.c \
test_read_vtk.c \
test_refine_by_size.c \
test_refine_topology.c \
test_up_from_down.c \
test_vtk.c \
test_warp.c \
test_warp_3d.c \
test_vtk_surfer.c \
test_vtkdiff.c \
test_push_cloud.c \
test_qr.c \
test_box.c \
test_pvtu.c \
test_split.c \
test_global_part.c \
test_loop.c \
test_base64.c

lib_sources := \
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
loop_host.c \
inertia.c \
derive_faces.c \
node_ele_io.c \
cloud.c \
tag.c \
form_cloud.c \
element_field.c \
mesh_diff.c \
push_cloud.c \
qr.c \
omega_h.c \
comm.c \
files.c \
global.c \
base64.c

#handle optional features:
USE_MPI ?= 0
#comm.c is compiled with -DUSE_MPI
objs/comm.o : CFLAGS += -DUSE_MPI=$(USE_MPI)

LOOP_MODE ?= serial
ifeq "$(LOOP_MODE)" "cuda"
lib_sources += loop_cuda.c
endif

#generated file names are derived from source
#file names by simple patterns:
exes := $(patsubst test_%.c,bin/%.exe,$(test_sources))
test_objects := $(patsubst %.c,objs/%.o,$(test_sources))
lib_objects := $(patsubst %.c,objs/%.o,$(lib_sources))
sources := $(lib_sources) $(test_sources)
depfiles := $(patsubst %.c,deps/%.dep,$(sources))

lib := lib/libomega_h.a

#the default compilation target is to compile
#the library and all executables
all: $(lib) $(exes)

#cleanup removes dependency files, object files,
#and executables
clean:
	rm -rf deps/*.dep objs/*.o bin/*.exe lib/*.a loop.h

#"all" and "clean" are targets, not files or directories
.PHONY: all clean

#our rule for compiling a source file to an
#object, specifies that the object goes in objs/
objs/%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

$(lib): $(lib_objects)
	ar cru $@ $(lib_objects)
	ranlib $@

#general rule for an executable: link its object
#file with all the $(common_objects)
# $@ is the thing being built and $^ is all
#the things it depends on (the objects)
bin/%.exe: objs/test_%.o $(lib)
	$(CC) $(LDFLAGS) -o $@ objs/test_$*.o $(lib) $(LDLIBS)

#loop.h is a copy of one of several existing files,
#chosen at compile time based on the kind of
#shared memory loop parallelism we want
loop.h : loop_$(LOOP_MODE).h
	cp $< $@

#Copied this mess from the GNU make documentation.
#It generates dependency files from source files,
#and uses SED to change the rules
#such that the output is both an object file and a
#dependency file.
#It warrants further explanation:
#  cc -MM foo.c
#will produce a dependency line such as:
#  foo.o : foo.c foo.h bar.h
#The SED script changes this to:
#  objs/foo.o deps/foo.dep : foo.c foo.h bar.h
#$* is the same as the % in the rule pattern,
#"foo" in the example above.
#the $@ will insert the "deps/foo.dep"

#loop.h is thrown in as a dependency because it
#may not exist when the depfiles are being generated,
#causing an error when cc -MM doesn't find it,
#and the knowedge that the depfile depends on it
#in the depfile itself !
deps/%.dep: %.c loop.h
	set -e; rm -f $@; \
	$(CPP) -M $(CPPFLAGS) $< > $@.in; \
	sed 's,$*\.o,objs/$*.o $@,g' < $@.in > $@; \
	rm -f $@.in

#include the auto-generated dependency files for
#all source files.
#this is the funny recursion that keeps
#header file dependencies worked out at all times.
include $(depfiles)
