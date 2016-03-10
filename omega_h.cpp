#include "include/omega_h.hpp"

#include "adapt.hpp"
#include "comm.hpp"
#include "ghost_mesh.hpp"
#include "loop.hpp"
#include "mark.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "size.hpp"
#include "tag.hpp"
#include "tables.hpp"
#include "vtk_io.hpp"

/*@
  osh_init - Initializes omega_h.

   If MPI was not already initialized, this
   will call MPI_Init.
   It currently does nothing else, but please
   consider it standard practice to call this
   before any other omega_h APIs.

  Collective

  Level: beginner

.seealso: osh_fini()
@*/
void osh_init(void)
{
  comm_init();
}

/*@
  osh_fini - Finalizes omega_h (opposite osh_init).

   If osh_init calls MPI_Init, this will call MPI_Finalize.
   It currently does nothing else, but please
   consider it standard practice to call this
   after all other omega_h APIs.

  Collective

  Level: beginner

.seealso: osh_init()
@*/
void osh_fini(void)
{
  comm_fini();
}

/*@
  osh_free - Destroys a mesh.

   All memory associated with the mesh
   is freed, and the mesh handle is invalid
   after this call.

  Collective

  Input Parameters:
. m - the mesh handle

  Level: beginner

.seealso: osh_new(), osh_read_vtk()
@*/
void osh_free(osh_t m)
{
  free_mesh((struct mesh*)m);
}

/*@
  osh_read_vtk - Read a mesh from our VTK format.

   This is the preferred way of creating a mesh
   data structure.
   See our conversion programs for converting
   meshes from Gmsh, TetGen, Triangle, etc. into
   VTK format.
   Also see our box program which generates a
   mesh of a simple unit domain.

   If a ".vtu" filename is given while running
   a single MPI rank, a serial mesh will be read
   from that file.
   If a ".vtu" filename is given while running more
   than one MPI rank, the serial mesh will subsequently
   be partitioned amongs the ranks.
   If a ".pvtu" filename is given, the pre-partitioned
   mesh will be loaded from the associated ".vtu" files,
   one file per MPI rank.

   The mesh new handle is returned.

  Collective

  Input Parameters:
. filename - the path to the mesh file(s), see description

  Level: beginner

.seealso: osh_write_vtk(), osh_free(), osh_new()
@*/
osh_t osh_read_vtk(char const* filename)
{
  return (osh_t) read_mesh_vtk(filename);
}

/*@
  osh_write_vtk - Write a mesh to our VTK format.

   If the mesh is parallel, a ".pvtu" extension
   should be used which creates a ".pvtu" file and
   corresponding ".vtu" files for each MPI rank.

   Otherwise, use the ".vtu" extension.

  Collective

  Input Parameters:
  + m - the mesh
  - filename - the output path, see description

  Level: beginner

.seealso: osh_read_vtk()
@*/
void osh_write_vtk(osh_t m, char const* filename)
{
  write_mesh_vtk((struct mesh*)m, filename);
}

/*@
  osh_new - Create a raw empty mesh.

   This is the advanced method of constructing
   meshes; it allows the user to define connectivity
   arrays directly.
   Don't use this unless you know about full
   topology representation.
   The mesh handle returned from this function is
   an empty mesh for which all dimensions
   have to be defined.

  Collective

  Input Parameters:
. elem_dim - Range [1,3], default 3, the dimension of elements

  Level: advanced

.seealso: osh_build_ents(), osh_dim()
@*/
osh_t osh_new(unsigned elem_dim)
{
  /* TODO: right now only APF uses this function,
     so its set to FULL/parallel for their convenience. */
  return (osh_t) new_mesh(elem_dim, MESH_FULL, 1);
}

/*@
  osh_build_ents - Create entities manually.

   This function allocates an array for the entity-to-vertex
   connectivity of the given dimension entities,
   sets it to represent connectivity in the mesh,
   and returns it to the user.

   The user should fill in the array right after
   receiving it from this function, and before
   calling any other API functions.
   See osh_down for the format of this array.

   The array is returned rather than received to
   allow omega_h to use its own allocators.

   In CUDA mode, the array will be in device memory.

  Collective

  Input Parameters:
+ m - mesh handle returned from osh_new
. ent_dim - dimension of the entities
- nents - the number of entities

  Level: advanced

.seealso: osh_new(), osh_down()
@*/
unsigned* osh_build_ents(osh_t m, unsigned ent_dim, unsigned nents)
{
  unsigned nverts_per_ent = 0;
  if (ent_dim)
    nverts_per_ent = the_down_degrees[ent_dim][0];
  unsigned* conn = LOOP_MALLOC(unsigned, nents * nverts_per_ent);
  /* this relies on the fact that mesh_set_ents doesn't expect
     the contents of "conn" to be filled in yet */
  mesh_set_ents((struct mesh*)m, ent_dim, nents, conn);
  return conn;
}

/*@
  osh_dim - Get the mesh / element dimension.

   Returns 1 for 1D meshes, 2 for 2D triangle meshes, etc.

  Input Parameters:
. m - mesh handle

  Level: beginner

.seealso: osh_count()
@*/
unsigned osh_dim(osh_t m)
{
  return mesh_dim((struct mesh*)m);
}

/*@
  osh_nelems - Get the number of mesh elements.

   Returns the number of elements on this MPI rank.
   Equivalent to osh_count(m, osh_dim(m)).

  Input Parameters:
. m - mesh handle

  Level: beginner

.seealso: osh_count()
@*/
unsigned osh_nelems(osh_t m)
{
  return mesh_count((struct mesh*)m,
      mesh_dim((struct mesh*)m));
}

/*@
  osh_nverts - Get the number of mesh vertices.

   Returns the number of vertices on this MPI rank.
   Equivalent to osh_count(m, 0).

  Input Parameters:
  . m - mesh handle

  Level: beginner

.seealso: osh_count()
@*/
unsigned osh_nverts(osh_t m)
{
  return mesh_count((struct mesh*)m, 0);
}

/*@
  osh_count - Get the number of mesh entities of one dimension

   Returns the number of entities of the given
   dimension on this MPI rank.

  Input Parameters:
. m - mesh handle

  Level: beginner

.seealso: osh_count()
@*/
unsigned osh_count(osh_t m, unsigned dim)
{
  return mesh_count((struct mesh*)m, dim);
}

/*@
  osh_down - Get the connectivity from high dimension to low

   Returns a read-only array describing connectivity
   from the high dimension to the lower one, on this MPI rank.

   The array is organized by higher entity,
   and for each higher entity there are contiguous
   entries pointing to lower entities.
   For example, from triangles to vertices, the array
   would be of size (3 * osh_count(m, 2)), and the entries

.vb
   unsigned const* down = osh_down(m, 2, 0);
   for (unsigned j = i * 3; j < (i + 1) * 3; ++j)
     down[j];
.ve

   contain the vertex indices of triangle (i).
   All indices are local to the MPI rank.

  Input Parameters:
+ m - mesh handle
. high_dim - the higher dimension
- low_dim - the lower dimension

  Level: beginner

.seealso: osh_count(), osh_dim(), osh_up(), osh_up_offs()
@*/
unsigned const* osh_down(osh_t m, unsigned high_dim, unsigned low_dim)
{
  return mesh_ask_down((struct mesh*)m, high_dim, low_dim);
}

/*@
  osh_up - Get the connectivity from low dimension to high (1)

   Returns a read-only array describing connectivity
   from the low dimension to the higher one, on this MPI rank.
   Works together with osh_up_offs().

   Similar to osh_down(), upward adjacencies are organized
   by lower entities, for each lower entity there
   are contiguous entries pointing to higher entities.
   However, the number of entries per lower entity is variable
   so there is a separate array of offsets describing where
   they start and stop.
   For example, from vertices to triangles, the array
   would be of size (osh_offs(m,0,2)[osh_count(m,0)]),
   and the entries

.vb
   unsigned const* offs = osh_up_offs(m,0,2);
   unsigned const* up = osh_up(m,0,2);
   for (unsigned j = offs[i]; j < offs[i + 1]; ++j)
     up[j];
.ve

   contain the vertex indices of triangle (i).
   All indices are local to the MPI rank.

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
. low_dim - the lower dimension
- high_dim - the higher dimension

  Level: beginner

.seealso: osh_up_offs(), osh_up_dir(), osh_down(), osh_dim(), osh_count()
@*/
unsigned const* osh_up(osh_t m, unsigned low_dim, unsigned high_dim)
{
  return mesh_ask_up((struct mesh*)m, low_dim, high_dim)->adj;
}

/*@
  osh_up_offs - Get the connectivity from low dimension to high (2)

   Returns a read-only array describing connectivity offsets
   from the low dimension to the higher one, on this MPI rank.
   Works together with osh_up().
   This array is of size (osh_count(m, low_dim) + 1).

   In CUDA mode, the array is in device memory.

  Input Parameters:
  + m - mesh handle
  . low_dim - the lower dimension
  - high_dim - the higher dimension

  Level: beginner

.seealso: osh_up(), osh_up_dir(), osh_down(), osh_dim(), osh_count()
@*/
unsigned const* osh_up_offs(osh_t m, unsigned low_dim, unsigned high_dim)
{
  return mesh_ask_up((struct mesh*)m, low_dim, high_dim)->offsets;
}

/*@
  osh_up_dirs - Get the downward order of upward connectivity

   Returns a read-only array that lines up with the osh_up() array.
   It indicates the local downward order for one connectivity relation
   ship.
   For example, say that triangle 42 has vertices {22,18,33}.
   Then, for the (18 -> 42) entry in osh_up(), this array contains (1).

   In other words:
.vb
   unsigned const* offs = osh_up_offs(m,0,2);
   unsigned const* up = osh_up(m,0,2);
   unsigned const* dir = osh_up_dir(m,0,2);
   unsigned const* down = osh_down(m,0,2);
   for (unsigned i = 0; i < osh_count(m,0); ++i)
     for (unsigned j = offs[i]; j < offs[i + 1]; ++j)
       assert(i == down[up[j] * 3 + dir[j]]);
.ve

  Input Parameters:
+ m - mesh handle
. low_dim - the lower dimension
- high_dim - the higher dimension

  Level: intermediate

.seealso: osh_up(), osh_up_offs(), osh_down(), osh_dim(), osh_count()
@*/
unsigned const* osh_up_dirs(osh_t m, unsigned low_dim, unsigned high_dim)
{
  return mesh_ask_up((struct mesh*)m, low_dim, high_dim)->directions;
}

/*@
  osh_star - Get the connectivity across other entities (1)

   Returns a read-only array describing connectivity
   from the low dimension across the higher one
   and back to the low dimension, on this MPI rank.
   Works together with osh_star_offs().
   The storage format is similar to osh_up() and osh_up_offs().

   The information is limited to this MPI rank.
   Specifically, if the mesh is not ghosted, entities
   at the partition boundary will have incomplete information.

   In CUDA mode, the array is in device memory.

  Examples are:

   osh_star(m, 0, 1) - vertices-to-vertices across edges

   osh_star(m, 0, 2) - vertices-to-vertices across triangles

   osh_star(m, 1, 3) - edges-to-edges across tetrahedra

  Input Parameters:
+ m - mesh handle
. low_dim - the lower dimension
- high_dim - the higher dimension

  Level: intermediate

.seealso: osh_star_offs(), osh_dim(), osh_count()
@*/
unsigned const* osh_star(osh_t m, unsigned low_dim, unsigned high_dim)
{
  return mesh_ask_star((struct mesh*)m, low_dim, high_dim)->adj;
}

/*@
  osh_star_offs - Get the connectivity across other entities (2)

   Returns a read-only array describing connectivity
   from the low dimension across the higher one
   and back to the low dimension, on this MPI rank.
   Works together with osh_star().

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
. low_dim - the lower dimension
- high_dim - the higher dimension

  Level: intermediate

.seealso: osh_star(), osh_dim(), osh_count()
@*/
unsigned const* osh_star_offs(osh_t m, unsigned low_dim, unsigned high_dim)
{
  return mesh_ask_star((struct mesh*)m, low_dim, high_dim)->offsets;
}

/*@
  osh_coords - Get the vertex coordinates.

   Returns a read-only array containing the {x,y,z}
   coordinates for each vertex, organized by vertex.
   This array has three entries per vertex regardless
   of mesh dimension.

   In CUDA mode, the array is in device memory.

  Input Parameters:
. m - mesh handle

  Level: beginner

.seealso: osh_down(), osh_count(), osh_nverts()
@*/
double const* osh_coords(osh_t m)
{
  return mesh_find_tag((struct mesh*)m, 0, "coordinates")->d.f64;
}

/*@
  osh_own_rank - Get the owning MPI ranks of entities.

   Returns a read-only array containing the owner
   MPI rank for each entity of a given dimension.
   This ownership is often determined by omega_h
   and affects the behavior of osh_accumulate_to_owner()
   and osh_conform().

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
- dim - the dimension of entities to inquire about

  Level: intermediate

.seealso: osh_own_id(), osh_accumulate_to_owner(), osh_conform(), osh_global()
@*/
unsigned const* osh_own_rank(osh_t m, unsigned dim)
{
  return mesh_ask_own_ranks((struct mesh*)m, dim);
}

/*@
  osh_own_id - Get the owning MPI ranks of entities.

   Returns a read-only array containing the owner
   local index for each entity of a given dimension.
   This ownership is often determined by omega_h
   and affects the behavior of osh_accumulate_to_owner()
   and osh_conform().
   The indices listed here are for the owner MPI
   rank, so they are usually only useful when transmitted
   to the owner rank.

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
- dim - the dimension of entities to inquire about

  Level: advanced

.seealso: osh_own_rank(), osh_accumulate_to_owner(), osh_conform(), osh_global()
@*/
unsigned const* osh_own_id(osh_t m, unsigned dim)
{
  return mesh_ask_own_ids((struct mesh*)m, dim);
}

/*@
  osh_global - Get the global numbers of entities.

   Returns a read-only array containing the global number
   for each entity of a given dimension.
   This numbering is often determined by omega_h.

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
- dim - the dimension of entities to inquire about

  Level: beginner

.seealso: osh_own_rank(), own_id()
@*/
unsigned long const* osh_global(osh_t m, unsigned dim)
{
  return mesh_ask_globals((struct mesh*)m, dim);
}

/*@
  osh_new_field - Create a floating-point field unless it exists.

   If a field by this name does not exist over entities
   of the given dimension, one will be created.
   Otherwise this is equivalent to osh_get_field().
   The internally allocated array of data is returned,
   and can be filled in after this call.
   The components for each entity are contiguous.

   The array is returned rather than received to
   allow omega_h to use its own allocators.

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
. dim - the dimension of entities
. name - the field name
- ncomps - the number of components (doubles) per entity

  Level: beginner

.seealso: osh_free_field(), osh_get_field()
@*/
double* osh_new_field(osh_t m, unsigned dim, char const* name,
    unsigned ncomps, enum osh_transfer tt)
{
  if (mesh_find_tag((struct mesh*)m, dim, name))
    return osh_get_field(m, dim, name);
  double* data = LOOP_MALLOC(double, ncomps * mesh_count((struct mesh*)m, dim));
  add_tag2(mesh_tags((struct mesh*)m, dim), TAG_F64, name, ncomps, tt, data);
  return data;
}

/*@
  osh_get_field - Access data for a field.

   Returns the internal data array of a field.

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
. dim - the dimension of entities
- name - the field name

  Level: beginner

.seealso: osh_new_field()
@*/
double* osh_get_field(osh_t m, unsigned dim, char const* name)
{
  return mesh_find_tag((struct mesh*)m, dim, name)->d.f64;
}

/*@
  osh_free_field - Destroy a field.

  Input Parameters:
+ m - mesh handle
. dim - the entity dimension over which the field is defined
- name - the field name

  Level: beginner

.seealso: osh_new_field()
@*/
void osh_free_field(osh_t m, unsigned dim, char const* name)
{
  mesh_free_tag((struct mesh*)m, dim, name);
}

/*@
  osh_nfields - Count the number of fields in a dimension.

   This is useful in conjunction with osh_field() to iterate
   over all fields.

  Input Parameters:
+ m - mesh handle
- dim - the dimension of entities

  Level: intermediate

.seealso: osh_field()
@*/
unsigned osh_nfields(osh_t om, unsigned dim)
{
  struct mesh* m = (struct mesh*) om;
  unsigned n = 0;
  for (unsigned i = 0; i < mesh_count_tags(m, dim); ++i)
    if (mesh_get_tag(m, dim, i)->type == TAG_F64)
      ++n;
  return n;
}

/*@
  osh_field - Access the i'th field.

   Returns the name of the i'th field in
   this dimension.
   The name can then be given to osh_get_field(), etc.

  Input Parameters:
+ m - mesh handle
. dim - the dimension of entities
- i - index of the field

  Level: intermediate

.seealso: osh_nfields(), osh_get_field()
@*/
char const* osh_field(osh_t om, unsigned dim, unsigned i)
{
  struct mesh* m = (struct mesh*) om;
  unsigned n = 0;
  for (unsigned j = 0; j < mesh_count_tags(m, dim); ++j)
    if (mesh_get_tag(m, dim, j)->type == TAG_F64)
      if (i == n++)
        return mesh_get_tag(m, dim, j)->name;
  return 0;
}

/*@
  osh_components - Get the number of components of a field or label.

  Input Parameters:
+ m - mesh handle
. dim - the dimension of entities
- name - the field or label name

  Level: intermediate

.seealso: osh_nfields()
@*/
unsigned osh_components(osh_t m, unsigned dim, char const* name)
{
  return mesh_find_tag((struct mesh*)m, dim, name)->ncomps;
}

/*@
  osh_new_label - Create an integer label unless it exists

   If a label with this name does not exist over entities
   of the given dimension, one will be created.
   Otherwise this is equivalent to osh_get_label().
   The internally allocated array of data is returned,
   and can be filled in after this call.
   The components for each entity are contiguous.

   The array is returned rather than received to
   allow omega_h to use its own allocators.

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
. dim - the dimension of entities
. name - the label name
- name - the number of components (integers)

  Level: beginner

.seealso: osh_get_label(), osh_free_label()
@*/
unsigned* osh_new_label(osh_t m, unsigned dim, char const* name, unsigned ncomps)
{
  if (mesh_find_tag((struct mesh*)m, dim, name))
    return osh_get_label(m, dim, name);
  unsigned* data = LOOP_MALLOC(unsigned,
      ncomps * mesh_count((struct mesh*)m, dim));
  mesh_add_tag((struct mesh*)m, dim, TAG_U32, name, ncomps, data);
  return data;
}

/*@
  osh_get_label - Access data for a label.

   Returns the internal data array of a label.

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
. dim - the dimension of entities
- name - the label name

  Level: beginner

.seealso: osh_new_label()
@*/
unsigned* osh_get_label(osh_t m, unsigned dim, char const* name)
{
  return mesh_find_tag((struct mesh*)m, dim, name)->d.u32;
}

/*@
  osh_free_label - Destroy a label.

  Input Parameters:
+ m - mesh handle
. dim - the entity dimension over which the label is defined
- name - the label name

  Level: beginner

.seealso: osh_new_label()
@*/
void osh_free_label(osh_t m, unsigned dim, char const* name)
{
  mesh_free_tag((struct mesh*)m, dim, name);
}

/*@
  osh_new_global - Create a global numbering of entities

   This should only be used in conjunction with
   osh_new() to create a mesh manually in parallel.
   The internally allocated array of global numbers is returned,
   and should be filled in after this call and before
   any subsequent API calls.

   The array is returned rather than received to
   allow omega_h to use its own allocators.

   In CUDA mode, the array is in device memory.

  Input Parameters:
+ m - mesh handle
- dim - the dimension of entities

  Level: advanced

.seealso: osh_global()
@*/
unsigned long* osh_new_global(osh_t m, unsigned dim)
{
  unsigned nents = mesh_count((struct mesh*)m, dim);
  unsigned long* global = LOOP_MALLOC(unsigned long, nents);
  mesh_set_globals((struct mesh*)m, dim, global);
  return global;
}

/*@
  osh_accumulate_to_owner - Sum field values between MPI ranks.

   For the specified field, omega_h will transmit
   value from non-owned entities to owned entities
   and add them to the owned entity values.
   Ownership is defined as in osh_own_rank() and osh_own_id().

  Input Parameters:
+ m - mesh handle
. dim - the entity dimension over which the field is defined
- name - the field name

  Level: intermediate

.seealso: osh_own_rank(), osh_own_id(), osh_conform()
@*/
void osh_accumulate_to_owner(osh_t m, unsigned dim, char const* name)
{
  mesh_accumulate_tag((struct mesh*)m, dim, name);
}

/*@
  osh_conform - Ensure consistent field values between MPI ranks.

   For the specified field, omega_h will transmit
   value from owned entities to non-owned entities
   and overwrite non-owned values.
   Ownership is defined as in osh_own_rank() and osh_own_id().

  Input Parameters:
+ m - mesh handle
. dim - the entity dimension over which the field is defined
- name - the field name

  Level: intermediate

.seealso: osh_own_rank(), osh_own_id(), osh_accumulate_to_owner()
@*/
void osh_conform(osh_t m, unsigned dim, char const* name)
{
  mesh_conform_tag((struct mesh*)m, dim, name);
}

LOOP_KERNEL(mark_or_kern,
    unsigned const* to_mark,
    unsigned* marked)
  if (to_mark[i])
    marked[i] = 1;
}

/*@
  osh_mark_classified - Mark entities classified on a certain boundary.

   Given a dimension and boundary classification specification
   in the form of boundary dimension and integer id,
   as well as a user-allocated array that indicates marked
   entities, this function will modify the marked array
   to include entities classified on the specified boundary.

   Boundary markers are either derived by omega_h or taken
   to be consistent with the Gmsh geometry if that is the
   origin of the mesh.

  Input Parameters:
+ m - mesh handle
. ent_dim - the dimension of the entities
. class_dim - the dimension of the boundary entity
. class_id - the ID of the boundary entity
- marked - user array, marked[i] = 0 or 1

  Level: intermediate
@*/
void osh_mark_classified(osh_t m, unsigned ent_dim,
    unsigned class_dim, unsigned class_id, unsigned* marked)
{
  unsigned* to_mark = mesh_mark_class((struct mesh*)m, ent_dim,
      class_dim, class_id);
  unsigned nents = osh_count(m, ent_dim);
  LOOP_EXEC(mark_or_kern, nents, to_mark, marked);
  loop_free(to_mark);
}

/*@
  osh_ghost - Create ghost/halo layers of elements.

   This function adds ghost copies of elements from other
   MPI ranks.
   Adding a layer begins by looking at all vertices
   (owned and not) which are visible on this MPI rank,
   and requesting all elements which use those vertices
   (and subsequently all entities on the boundary of those elements).
   This is done recursively until the desired number
   of layers is established.
   Ghost elements may be identified by osh_own_rank(),
   they will have owner MPI ranks other than the local one.

   If the desired number of layers already exists,
   no change occurs.
   One may specify nlayers=0 to remove ghost layers.

   This operation modifies all connectivity / field
   arrays, so any user information that relies on
   their pointers or content must be reconstructed
   afterwards.

  Input Parameters:
+ m - mesh handle
- nlayers - number of ghost layers to add

  Level: advanced
@*/
void osh_ghost(osh_t m, unsigned nlayers)
{
  ghost_mesh((struct mesh*)m, nlayers);
}

/*@
  osh_adapt - Adapts the mesh to a given size field.

   The size field is currently a single value per vertex
   describing the desired length of mesh edges at that location,
   it must be called "adapt_size".
   The return value will be zero if no connectivity changes are made.

  Collective

  Input Parameters:
+ m - The mesh object.
. size_quality_floor - Range [0.0 - 0.5], default 1/3. Edges shorter than this
                       times their desired size will attempt coarsening.
. good_element_quality - Range [0.0 - 1.0], default 0.3. No elements will lower
                         quality will be created, and existing elements with
                         lower quality will attempt shape correction.
. nsliver_layers - Range [0 - ?], default 4. Number of layers of elements around
                   a sliver to consider for correction.
- max_passes - Range [1 - 100], default 50. Maximum number of modification passes
               before the program fails.

  Level: advanced

.keywords: adapt
.seealso: osh_new_field(), osh_identity_size()
@*/
unsigned osh_adapt(osh_t m,
    double size_ratio_floor,
    double good_element_quality,
    unsigned nsliver_layers,
    unsigned max_passes)
{
  return mesh_adapt((struct mesh*)m,
      size_ratio_floor,
      good_element_quality,
      nsliver_layers,
      max_passes);
}

/*@
  osh_identity_size - Compute the "current" size field.

  This function computes a size field which,
  when used by osh_adapt(), will cause little to no
  change in the mesh edge lengths.
  To be precise, it computes at each vertex the
  length of the longest adjacent edge.

  There are several reasons why this might be useful.
  First, one can use this to limit adaptation to
  just shape correction, and improve the quality of
  elements while maintaining resolution.

  One can also compute this field before Langrangian
  mesh motion (moving vertices) and later use osh_adapt
  to restore the previous resolution while roughly
  maintaining new vertex locations.

  Finally, the identity size field may be used as the
  starting point for the final size field, for example
  it may be divided by some factor in order to increase
  resolution throughout by that factor, even if the
  initial resolution was not uniform.

  Input Parameters:
+ m - mesh handle
- name - the name used for the resulting field

  Level: advanced
@*/
void osh_identity_size(osh_t m, char const* name)
{
  mesh_identity_size_field((struct mesh*)m, name);
}
