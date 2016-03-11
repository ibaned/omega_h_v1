#The user configuration is placed in config.mk
#Users edit config.mk but should never need
#to edit the Makefile.
#This also separates the build configuration
#in the eye of the version control system,
#which does not keep track of config.mk
include config.mk

test_sources := \
test_qr.cpp \
test_reorder.cpp \
test_identity.cpp \
test_flounder.cpp \
test_one_refine.cpp \
test_one_coarsen.cpp \
test_one_swap.cpp \
test_partition.cpp \
test_box.cpp \
test_node_ele.cpp \
test_node_ele_attrib.cpp \
test_from_gmsh.cpp \
test_vtk_ascii.cpp \
test_vtkdiff.cpp \
test_coarsen_by_size.cpp \
test_refine_by_size.cpp \
test_grad.cpp \
test_warp.cpp \
test_warp_3d.cpp \
test_warp_perf.cpp \
test_migrate.cpp \
test_conform.cpp \
test_ghost.cpp \
test_memory.cpp \
test_subdim.cpp \
test_loop.cpp \
test_to_la.cpp

lib_sources := \
swap_fit.cpp \
coarsen_fit.cpp \
refine_fit.cpp \
shuffle_mesh.cpp \
reorder.cpp \
bfs.cpp \
star.cpp \
tables.cpp \
up_from_down.cpp \
ints.cpp \
vtk_io.cpp \
refine_topology.cpp \
splits_to_domains.cpp \
quality.cpp \
size.cpp \
bridge_graph.cpp \
refine_common.cpp \
refine.cpp \
indset.cpp \
reflect_down.cpp \
dual.cpp \
refine_nodal.cpp \
refine_conserve.cpp \
refine_qualities.cpp \
doubles.cpp \
coarsen.cpp \
check_collapse_class.cpp \
coarsen_qualities.cpp \
coarsen_topology.cpp \
coarsen_conserve.cpp \
collapses_to_verts.cpp \
collapses_to_ents.cpp \
mesh.cpp \
graph.cpp \
warp_to_limit.cpp \
eval_field.cpp \
element_gradients.cpp \
recover_by_volume.cpp \
size_from_hessian.cpp \
subset.cpp \
adapt.cpp \
coarsen_common.cpp \
mark.cpp \
swap.cpp \
swap_qualities.cpp \
swap_topology.cpp \
swap_conserve.cpp \
edge_ring.cpp \
edge_swap.cpp \
edge_swap_edges.cpp \
loop_host.cpp \
inertia.cpp \
derive_sides.cpp \
node_ele_io.cpp \
tag.cpp \
element_field.cpp \
mesh_diff.cpp \
qr.cpp \
omega_h.cpp \
comm.cpp \
files.cpp \
global.cpp \
base64.cpp \
invert_map.cpp \
owners_from_global.cpp \
gmsh_io.cpp \
exchanger.cpp \
parallel_inertial_bisect.cpp \
parallel_mesh.cpp \
parallel_modify.cpp \
arrays.cpp \
migrate_mesh.cpp \
bcast.cpp \
close_partition.cpp \
ghost_mesh.cpp \
derive_model.cpp \
compress.cpp \
inherit.cpp

#handle optional features:
PREFIX ?= /usr/local
USE_ZLIB ?= 0
USE_MPI ?= 0
USE_MPI3 ?= $(USE_MPI)
USE_CUDA_MALLOC_MANAGED ?= 0
MEASURE_MEMORY ?= 0
LOOP_MODE ?= serial
MPIRUN ?= mpirun
VALGRIND ?= ""
PATIENT ?= 0
SHARED ?= 0
#comm.cpp is compiled with -DUSE_MPI=
objs/comm.o : CPPFLAGS += -DUSE_MPI=$(USE_MPI)
deps/comm.dep : CPPFLAGS += -DUSE_MPI=$(USE_MPI)
ifeq "$(USE_MPI)" "1"
  lib_sources += compat_mpi.cpp
  objs/compat_mpi.o : CPPFLAGS += -DUSE_MPI3=$(USE_MPI3)
  deps/compat_mpi.dep : CPPFLAGS += -DUSE_MPI3=$(USE_MPI3)
endif
objs/loop_host.o : CPPFLAGS += -DMEASURE_MEMORY=$(MEASURE_MEMORY)
ifeq "$(LOOP_MODE)" "cuda"
  lib_sources += loop_cuda.cpp
  objs/loop_cuda.o : CPPFLAGS += -DUSE_CUDA_MALLOC_MANAGED=$(USE_CUDA_MALLOC_MANAGED)
else
  test_sources += test_print_swap_edges.cpp
endif
ifeq "$(LOOP_MODE)" "kokkos"
  lib_sources += loop_kokkos.cpp
endif
objs/compress.o : CPPFLAGS += -DUSE_ZLIB=$(USE_ZLIB)
ifeq "$(USE_ZLIB)" "1"
  objs/compress.o : CPPFLAGS += -I$(ZLIB_INCLUDE)
endif
libraries := lib/libomega_h.a
ifeq "$(SHARED)" "1"
  CXXFLAGS += -fPIC -fvisibility=hidden
  libraries += lib/libomega_h.so
endif

#generated file names are derived from source
#file names by simple patterns:
exes := $(patsubst test_%.cpp,bin/%.exe,$(test_sources))
test_objects := $(patsubst %.cpp,objs/%.o,$(test_sources))
lib_objects := $(patsubst %.cpp,objs/%.o,$(lib_sources))
depfiles := $(patsubst %.cpp,deps/%.dep,$(lib_sources)) \
$(patsubst %.cpp,deps/%.dep,$(test_sources))

#the default compilation target is to compile
#the library and all executables
all: $(libraries) $(exes)

#cleanup removes dependency files, object files,
#and executables
clean:
	rm -rf deps/ objs/ bin/ lib/ loop.hpp

#just targets, not files or directories
.PHONY: all clean check install dep coverage

#our rule for compiling a source file to an
#object, specifies that the object goes in objs/
objs/%.o: %.cpp | objs
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

lib/libomega_h_internal.a: $(lib_objects) | lib
	ar cru $@ $(lib_objects)
	ranlib $@

lib/libomega_h.a: lib/libomega_h_internal.a
	cp $< $@

lib/libomega_h.so: $(lib_objects) | lib
	$(CXX) $(LDFLAGS) -shared -fPIC -o $@ $(lib_objects) $(LDLIBS)

#general rule for an executable: link its object
#file with the library
bin/%.exe: objs/test_%.o lib/libomega_h_internal.a | bin
	$(CXX) $(LDFLAGS) -L./lib -o $@ objs/test_$*.o -lomega_h_internal $(LDLIBS)

#loop.hpp is a copy of one of several existing files,
#chosen at compile time based on the kind of
#shared memory loop parallelism we want
loop.hpp : loop_$(LOOP_MODE).hpp
	cp $< $@

#make output directories if they don't yet exist.
#rules that generate files in these directories
#will depend on the directory as an "order-only"
#target by putting it after a pipe like this:
#  objs/%.o : %.cpp | objs
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
#  cc -M foo.cpp
#will produce a dependency line such as:
#  foo.o : foo.cpp foo.hpp bar.hpp
#The SED script changes this to:
#  objs/foo.o deps/foo.dep : foo.cpp foo.hpp bar.hpp
#$* is the same as the % in the rule pattern,
#"foo" in the example above.
#the $@ will insert the "deps/foo.dep"
#
#loop.hpp is thrown in as a dependency because it
#may not exist when the depfiles are being generated,
#causing an error when cc -M doesn't find it,
#and the knowedge that the depfile depends on it
#in the depfile itself !
deps/%.dep: %.cpp loop.hpp | deps
	set -e; rm -f $@; \
	$(CXX) -M $(CPPFLAGS) $< > $@.in; \
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
	install -m 644 $(libraries) $(PREFIX)/lib
	install -d $(PREFIX)/include
	install -m 644 include/omega_h.hpp $(PREFIX)/include
	install -d $(PREFIX)/share/man/man3
	install -m 644 share/man/man3/*.3 $(PREFIX)/share/man/man3

check: $(exes) data gold scratch
	MPIRUN=$(MPIRUN) VALGRIND=$(VALGRIND) \
  USE_MPI=$(USE_MPI) PATIENT=$(PATIENT) \
  LOOP_MODE=$(LOOP_MODE) ./run_tests.sh

data:
	git clone https://github.com/ibaned/omega_h_data.git data
gold:
	mkdir gold
scratch:
	mkdir scratch

coverage: objs scratch
	lcov --capture --directory objs --output-file scratch/coverage.info
	genhtml scratch/coverage.info --output-directory lcov-output

doc:
	doctext -mpath share/man/man3 -ext 3 -nolocation omega_h.cpp
