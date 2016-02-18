#The user configuration is placed in config.mk
#Users edit config.mk but should never need
#to edit the Makefile.
#This also separates the build configuration
#in the eye of the version control system,
#which does not keep track of config.mk
include config.mk

test_sources := \
test_identity.c \
test_flounder.c \
test_one_refine.c \
test_one_coarsen.c \
test_one_swap.c \
test_partition.c \
test_box.c \
test_node_ele.c \
test_from_gmsh.c \
test_vtk_ascii.c \
test_vtk_surfer.c \
test_vtkdiff.c \
test_coarsen_by_size.c \
test_refine_by_size.c \
test_grad.c \
test_warp.c \
test_warp_3d.c \
test_migrate.c \
test_conform.c \
test_ghost.c \
test_form_cloud.c \
test_fusion_part.c \
test_memory.c \
test_ask_up.c \
test_ask_down.c \
test_derive_model.c \
test_subdim.c \
test_loop.c

lib_sources := \
star.c \
tables.c \
up_from_down.c \
ints.c \
vtk.c \
refine_topology.c \
splits_to_domains.c \
quality.c \
size.c \
bridge_graph.c \
refine_common.c \
refine.c \
indset.c \
reflect_down.c \
dual.c \
refine_nodal.c \
refine_conserve.c \
refine_qualities.c \
doubles.c \
coarsen.c \
check_collapse_class.c \
coarsen_qualities.c \
coarsen_topology.c \
coarsen_conserve.c \
collapses_to_verts.c \
collapses_to_ents.c \
mesh.c \
graph.c \
warp_to_limit.c \
eval_field.c \
element_gradients.c \
recover_by_volume.c \
size_from_hessian.c \
subset.c \
adapt.c \
coarsen_common.c \
mark.c \
swap.c \
swap_qualities.c \
swap_topology.c \
swap_conserve.c \
edge_ring.c \
edge_swap.c \
edge_swap_edges.c \
loop_host.c \
inertia.c \
derive_sides.c \
node_ele_io.c \
cloud.c \
tag.c \
form_cloud.c \
element_field.c \
mesh_diff.c \
qr.c \
omega_h.c \
comm.c \
files.c \
global.c \
base64.c \
invert_map.c \
owners_from_global.c \
gmsh_io.c \
exchanger.c \
copy_tags.c \
parallel_inertial_bisect.c \
parallel_mesh.c \
parallel_modify.c \
migrate_cloud.c \
arrays.c \
migrate_mesh.c \
bcast.c \
close_partition.c \
ghost_mesh.c \
derive_model.c \
compress.c \
inherit.c

#handle optional features:
USE_ZLIB ?= 0
USE_MPI ?= 0
USE_MPI3 ?= 0
USE_CUDA_MALLOC_MANAGED ?= 0
MEASURE_MEMORY ?= 0
LOOP_MODE ?= serial
MPIRUN ?= mpirun
VALGRIND ?= ""
PATIENT ?= 0
#comm.c is compiled with -DUSE_MPI=
objs/comm.o : CPPFLAGS += -DUSE_MPI=$(USE_MPI)
deps/comm.dep : CPPFLAGS += -DUSE_MPI=$(USE_MPI)
ifeq "$(USE_MPI)" "1"
lib_sources += compat_mpi.c
objs/compat_mpi.o : CPPFLAGS += -DUSE_MPI3=$(USE_MPI3)
deps/compat_mpi.dep : CPPFLAGS += -DUSE_MPI3=$(USE_MPI3)
endif
objs/loop_host.o : CPPFLAGS += -DMEASURE_MEMORY=$(MEASURE_MEMORY)
lib_sources += loop_$(LOOP_MODE).c
ifeq "$(LOOP_MODE)" "cuda"
objs/loop_cuda.o : CPPFLAGS += -DUSE_CUDA_MALLOC_MANAGED=$(USE_CUDA_MALLOC_MANAGED)
else
test_sources += test_print_swap_edges.c
endif
objs/compress.o : CPPFLAGS += -DUSE_ZLIB=$(USE_ZLIB)
ifeq "$(USE_ZLIB)" "1"
objs/compress.o : CPPFLAGS += -I$(ZLIB_INCLUDE)
endif

#generated file names are derived from source
#file names by simple patterns:
exes := $(patsubst test_%.c,bin/%.exe,$(test_sources))
test_objects := $(patsubst %.c,objs/%.o,$(test_sources))
lib_objects := $(patsubst %.c,objs/%.o,$(lib_sources))
depfiles := $(patsubst %.c,deps/%.dep,$(lib_sources)) \
$(patsubst %.c,deps/%.dep,$(test_sources))

lib := lib/libomega_h.a

#the default compilation target is to compile
#the library and all executables
all: $(lib) $(exes)

#cleanup removes dependency files, object files,
#and executables
clean:
	rm -rf deps/ objs/ bin/ lib/ loop.h

#just targets, not files or directories
.PHONY: all clean check install dep

#our rule for compiling a source file to an
#object, specifies that the object goes in objs/
objs/%.o: %.c | objs
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

$(lib): $(lib_objects) | lib
	ar cru $@ $(lib_objects)
	ranlib $@

#general rule for an executable: link its object
#file with all the $(common_objects)
# $@ is the thing being built and $^ is all
#the things it depends on (the objects)
bin/%.exe: objs/test_%.o $(lib) | bin
	$(CC) $(LDFLAGS) -o $@ objs/test_$*.o $(lib) $(LDLIBS)

#loop.h is a copy of one of several existing files,
#chosen at compile time based on the kind of
#shared memory loop parallelism we want
loop.h : loop_$(LOOP_MODE).h
	cp $< $@

#make output directories if they don't yet exist.
#rules that generate files in these directories
#will depend on the directory as an "order-only"
#target by putting it after a pipe like this:
#  objs/%.o : %.c | objs
#this is needed because updating a file in a
#directory makes that directory "newer" than the
#file, which would cause an infinite loop in make.
objs:
	mkdir objs
bin:
	mkdir bin
lib:
	mkdir lib
deps:
	mkdir deps

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
#
#loop.h is thrown in as a dependency because it
#may not exist when the depfiles are being generated,
#causing an error when cc -MM doesn't find it,
#and the knowedge that the depfile depends on it
#in the depfile itself !
deps/%.dep: %.c loop.h | deps
	set -e; rm -f $@; \
	$(CPP) -M $(CPPFLAGS) $< > $@.in; \
	sed 's,$*\.o,objs/$*.o $@,g' < $@.in > $@; \
	rm -f $@.in

#a make target to just generate the dependency
#files without compiling source code yet
dep: $(depfiles)

#include the auto-generated dependency files for
#all source files.
#this is the funny recursion that keeps
#header file dependencies worked out at all times.
#the minus sign silences warnings when the
#depfiles don't exist yet.
-include $(depfiles)

install: all
	install -d $(PREFIX)/bin
	install -m 755 $(exes) $(PREFIX)/bin
	install -d $(PREFIX)/lib
	install -m 644 $(lib) $(PREFIX)/lib
	install -d $(PREFIX)/include
	install -m 644 include/omega_h.h $(PREFIX)/include

check: data $(exes)
	MPIRUN=$(MPIRUN) VALGRIND=$(VALGRIND) \
  USE_MPI=$(USE_MPI) PATIENT=$(PATIENT) \
  LOOP_MODE=$(LOOP_MODE) ./run_tests.sh

data:
	git clone https://github.com/ibaned/omega_h_data.git data
