#include "vtk.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arrays.h"
#include "base64.h"
#include "cloud.h"
#include "comm.h"
#include "files.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"
#include "tables.h"
#include "tag.h"

enum cell_type {
  VTK_VERTEX         = 1,
  VTK_POLY_VERTEX    = 2,
  VTK_LINE           = 3,
  VTK_POLY_LINE      = 4,
  VTK_TRIANGLE       = 5,
  VTK_TRIANGLE_STRIP = 6,
  VTK_POLYGON        = 7,
  VTK_PIXEL          = 8,
  VTK_QUAD           = 9,
  VTK_TETRA          =10,
  VTK_VOXEL          =11,
  VTK_HEXAHEDRON     =12,
  VTK_WEDGE          =13,
  VTK_PYRAMID        =14,
};

static enum cell_type const simplex_types[4] = {
  VTK_VERTEX,
  VTK_LINE,
  VTK_TRIANGLE,
  VTK_TETRA
};

static char const* type_name(enum tag_type t)
{
  switch (t) {
    case TAG_U8:  return "UInt8";
    case TAG_U32: return "UInt32";
    case TAG_U64: return "UInt64";
    case TAG_F64: return "Float64";
    default: return "";
  }
}

static char const* format_name(enum vtk_format fmt)
{
  switch (fmt) {
    case VTK_ASCII: return "ascii";
    case VTK_BINARY: return "binary";
    default: return "";
  }
}

static unsigned try_read_attrib(char const* elem, char const* name,
    char* val)
{
  char const* pname = strstr(elem, name);
  if (!pname)
    return 0;
  line_t assign;
  assert(strlen(pname) < sizeof(assign));
  strcpy(assign, pname);
  assert(assign[strlen(name) + 1] == '\"');
  char const* pval = strtok(assign + strlen(name) + 2, "\"");
  assert(pval && strlen(pval));
  strcpy(val, pval);
  return 1;
}

static void read_attrib(char const* elem, char const* name,
    char* val)
{
  unsigned ok = try_read_attrib(elem, name, val);
  assert(ok);
}

static void read_array_name(char const* header, char* name)
{
  read_attrib(header, "Name", name);
}

static enum tag_type read_array_type(char const* header)
{
  line_t val;
  read_attrib(header, "type", val);
  for (unsigned type = 0; type < TAG_TYPES; ++type)
    if (!strcmp(type_name((enum tag_type) type), val))
      return (enum tag_type) type;
  assert(0);
#ifdef __CUDACC__
  return TAG_U32;
#endif
}

static enum vtk_format read_array_format(char const* header)
{
  line_t val;
  read_attrib(header, "format", val);
  for (unsigned fmt = 0; fmt < VTK_FORMATS; ++fmt)
    if (!strcmp(format_name((enum vtk_format) fmt), val))
      return (enum vtk_format) fmt;
  assert(0);
#ifdef __CUDACC__
  return VTK_ASCII;
#endif
}

static unsigned read_int_attrib(char const* header, char const* attrib)
{
  line_t val;
  read_attrib(header, attrib, val);
  return (unsigned) atoi(val);
}

static unsigned try_read_int_attrib(char const* header, char const* attrib,
    unsigned* val)
{
  line_t val_text;
  unsigned ok = try_read_attrib(header, attrib, val_text);
  if (!ok)
    return 0;
  *val = (unsigned) atoi(val_text);
  return 1;
}

static unsigned read_array_ncomps(char const* header)
{
  unsigned n;
  if (try_read_int_attrib(header, "NumberOfComponents", &n))
    return n;
  return 1;
}

static void write_binary_array(FILE* file, enum tag_type t, unsigned nents,
    unsigned ncomps, void const* data)
{
  unsigned tsize = tag_size(t);
  unsigned size = tsize * ncomps * nents;
  char* s = base64_encode(&size, sizeof(size));
  fputs(s, file);
  loop_host_free(s);
  s = base64_encode(data, size);
  fputs(s, file);
  loop_host_free(s);
  fputc('\n', file);
}

static void* read_binary_array(FILE* file, enum tag_type t, unsigned nents,
    unsigned ncomps)
{
  unsigned long enc_nchars;
  char* enc = base64_fread(file, &enc_nchars);
  char const* p = enc;
  void* dec = base64_decode(&p, sizeof(unsigned));
  loop_host_free(dec);
  unsigned tsize = tag_size(t);
  dec = base64_decode(&p, nents * ncomps * tsize);
  loop_host_free(enc);
  void* swapped = generic_swap_if_needed(endianness(), nents * ncomps, tsize,
      dec);
  loop_host_free(dec);
  return swapped;
}

static void write_ascii_array(FILE* file, enum tag_type t, unsigned nents,
    unsigned ncomps, void const* data)
{
  switch (t) {
    case TAG_U8: {
      unsigned char const* p = (unsigned char const*) data;
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < ncomps; ++j)
          fprintf(file, " %hhu", *p++);
        fprintf(file, "\n");
      }
      break;
    }
    case TAG_U32: {
      unsigned const* p = (unsigned const*) data;
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < ncomps; ++j)
          fprintf(file, " %u", *p++);
        fprintf(file, "\n");
      }
      break;
    }
    case TAG_U64: {
      unsigned long const* p = (unsigned long const*) data;
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < ncomps; ++j)
          fprintf(file, " %lu", *p++);
        fprintf(file, "\n");
      }
      break;
    }
    case TAG_F64: {
      double const* p = (double const*) data;
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < ncomps; ++j)
          fprintf(file, " %.15e", *p++);
        fprintf(file, "\n");
      }
      break;
    }
  }
}

static void* read_ascii_array(FILE* file, enum tag_type type, unsigned nents,
    unsigned ncomps)
{
  unsigned n = nents * ncomps;
  switch (type) {
    case TAG_U8: {
      unsigned char* out = LOOP_HOST_MALLOC(unsigned char, n);
      for (unsigned i = 0; i < n; ++i)
        safe_scanf(file, 1, "%hhu", &out[i]);
      return out;
    }
    case TAG_U32: {
      unsigned* out = LOOP_HOST_MALLOC(unsigned, n);
      for (unsigned i = 0; i < n; ++i)
        safe_scanf(file, 1, "%u", &out[i]);
      return out;
    }
    case TAG_U64: {
      unsigned long* out = LOOP_HOST_MALLOC(unsigned long, n);
      for (unsigned i = 0; i < n; ++i)
        safe_scanf(file, 1, "%lu", &out[i]);
      return out;
    }
    case TAG_F64: {
      double* out = LOOP_HOST_MALLOC(double, n);
      for (unsigned i = 0; i < n; ++i)
        safe_scanf(file, 1, "%lf", &out[i]);
      return out;
    }
    default: return 0;
  }
}

static void describe_array(FILE* file, enum tag_type t,
    char const* name, unsigned ncomps, enum vtk_format fmt)
{
  fprintf(file, "type=\"%s\" Name=\"%s\""
      " NumberOfComponents=\"%u\" format=\"%s\"",
      type_name(t), name, ncomps, format_name(fmt));
}

static void describe_tag(FILE* file, struct const_tag* tag)
{
  describe_array(file, tag->type, tag->name, tag->ncomps, VTK_ASCII);
}

static void write_array(FILE* file, enum tag_type t,
    char const* name, unsigned nents, unsigned ncomps,
    void const* data, enum vtk_format fmt)
{
  fprintf(file, "<DataArray ");
  describe_array(file, t, name, ncomps, fmt);
  fprintf(file, ">\n");
  switch (fmt) {
    case VTK_ASCII:
      write_ascii_array(file, t, nents, ncomps, data);
      break;
    case VTK_BINARY:
      write_binary_array(file, t, nents, ncomps, data);
      break;
  }
  fprintf(file, "</DataArray>\n");
}

static unsigned seek_prefix_next(FILE* f,
    char* line, unsigned line_size, char const* prefix)
{
  unsigned long pl = strlen(prefix);
  if (!fgets(line, (int) line_size, f))
    return 0;
  return !strncmp(line, prefix, pl);
}

static void read_array(FILE* f, char const* line, enum tag_type* type,
    char* name, unsigned nents, unsigned* ncomps,
    void** data)
{
  *type = read_array_type(line);
  read_array_name(line, name);
  *ncomps = read_array_ncomps(line);
  enum vtk_format fmt = read_array_format(line);
  switch (fmt) {
    case VTK_ASCII:
      *data = read_ascii_array(f, *type, nents, *ncomps);
      break;
    case VTK_BINARY:
      *data = read_binary_array(f, *type, nents, *ncomps);
      break;
  }
  line_t tmpline;
  seek_prefix(f, tmpline, sizeof(tmpline), "</DataArray");
}

static unsigned read_tag(FILE* f, struct tags* ts, unsigned n)
{
  line_t line;
  if (!seek_prefix_next(f, line, sizeof(line), "<DataArray"))
    return 0;
  enum tag_type type;
  line_t name;
  unsigned ncomps;
  void* data;
  read_array(f, line, &type, name, n, &ncomps, &data);
  add_tag(ts, type, name, ncomps, data);
  return 1;
}

static void write_tag(FILE* file, unsigned nents, struct const_tag* tag,
    enum vtk_format fmt)
{
  write_array(file, tag->type, tag->name, nents, tag->ncomps, tag->d.raw,
      fmt);
}

static void write_connectivity(FILE* file, struct mesh* m, unsigned dim, enum vtk_format fmt)
{
  unsigned nents = mesh_count(m, dim);
  unsigned nverts_per_ent = the_down_degrees[dim][0];
  unsigned const* verts_of_elems = mesh_ask_down(m, dim, 0);
  write_array(file, TAG_U32, "connectivity", nents * nverts_per_ent, 1, verts_of_elems,
      fmt);
}

static void write_cell_arrays(FILE* file, struct mesh* m, enum vtk_format fmt)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  write_connectivity(file, m, elem_dim, fmt);
  unsigned* off = LOOP_HOST_MALLOC(unsigned, nelems);
  unsigned nverts_per_elem = the_down_degrees[elem_dim][0];
  for (unsigned i = 0; i < nelems; ++i)
    off[i] = (i + 1) * nverts_per_elem;
  write_array(file, TAG_U32, "offsets", nelems, 1, off, fmt);
  loop_host_free(off);
  unsigned char* types = LOOP_HOST_MALLOC(unsigned char, nelems);
  unsigned char type = (unsigned char) simplex_types[elem_dim];
  for (unsigned i = 0; i < nelems; ++i)
    types[i] = type;
  write_array(file, TAG_U8, "types", nelems, 1, types, fmt);
  loop_host_free(types);
}

static void write_tags(FILE* file, struct mesh* m, unsigned dim,
    enum vtk_format fmt, struct const_tag* except)
{
  for (unsigned i = 0; i < mesh_count_tags(m, dim); ++i) {
    struct const_tag* tag = mesh_get_tag(m, dim, i);
    if (tag != except)
      write_tag(file, mesh_count(m, dim), tag, fmt);
  }
}

static void write_unstructured_header(FILE* file)
{
  fprintf(file, "<VTKFile type=\"UnstructuredGrid\" byte_order=\"");
  if (endianness() == MY_LITTLE_ENDIAN)
    fprintf(file, "LittleEndian");
  else
    fprintf(file, "BigEndian");
  fprintf(file, "\">\n");
}

void write_vtu_opts(struct mesh* m, char const* filename, enum vtk_format fmt)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nverts = mesh_count(m, 0);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned do_edges = ((elem_dim > 1) && mesh_has_dim(m, 1));
  unsigned do_faces = ((elem_dim > 2) && mesh_has_dim(m, 2));
  FILE* file = fopen(filename, "w");
  write_unstructured_header(file);
  fprintf(file, "<UnstructuredGrid>\n");
  fprintf(file, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\"", nverts, nelems);
  if (do_edges)
    fprintf(file, " NumberOfEdges=\"%u\"", mesh_count(m, 1));
  if (do_faces)
    fprintf(file, " NumberOfFaces=\"%u\"", mesh_count(m, 2));
  fprintf(file, ">\n");
  fprintf(file, "<Points>\n");
  struct const_tag* coord_tag = mesh_find_tag(m, 0, "coordinates");
  write_tag(file, nverts, coord_tag, fmt);
  fprintf(file, "</Points>\n");
  if (do_edges) {
    fprintf(file, "<Edges>\n");
    write_connectivity(file, m, 1, fmt);
    fprintf(file, "</Edges>\n");
  }
  if (do_faces) {
    fprintf(file, "<Faces>\n");
    write_connectivity(file, m, 2, fmt);
    fprintf(file, "</Faces>\n");
  }
  fprintf(file, "<Cells>\n");
  write_cell_arrays(file, m, fmt);
  fprintf(file, "</Cells>\n");
  fprintf(file, "<PointData>\n");
  write_tags(file, m, 0, fmt, coord_tag);
  fprintf(file, "</PointData>\n");
  if (do_edges) {
    fprintf(file, "<EdgeData>\n");
    write_tags(file, m, 1, fmt, 0);
    fprintf(file, "</EdgeData>\n");
  }
  if (do_faces) {
    fprintf(file, "<FaceData>\n");
    write_tags(file, m, 2, fmt, 0);
    fprintf(file, "</FaceData>\n");
  }
  fprintf(file, "<CellData>\n");
  write_tags(file, m, elem_dim, fmt, 0);
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n");
  fprintf(file, "</UnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");
  fclose(file);
}

void write_vtu(struct mesh* m, char const* filename)
{
  write_vtu_opts(m, filename, VTK_BINARY);
}

static char const* the_step_prefix = 0;
static unsigned the_step = 0;

void start_vtk_steps(char const* prefix)
{
  the_step_prefix = prefix;
  the_step = 0;
}

void write_vtk_step(struct mesh* m)
{
  char fname[64];
  sprintf(fname, "%s_%04u.vtu", the_step_prefix, the_step);
  write_vtu(m, fname);
  ++the_step;
}

static void read_ent_counts(FILE* f, unsigned* nverts, unsigned* nelems,
    unsigned* do_edges, unsigned* do_faces,
    unsigned* nedges, unsigned* nfaces)
{
  line_t line;
  seek_prefix(f, line, sizeof(line), "<Piece");
  *nverts = read_int_attrib(line, "NumberOfPoints");
  *nelems = read_int_attrib(line, "NumberOfCells");
  *do_edges = try_read_int_attrib(line, "NumberOfEdges", nedges);
  *do_faces = try_read_int_attrib(line, "NumberOfFaces", nfaces);
}

static void read_nverts(FILE* f, unsigned* nverts)
{
  unsigned ignore;
  read_ent_counts(f, nverts, &ignore, &ignore, &ignore, &ignore, &ignore);
}

static unsigned read_dimension(FILE* f, unsigned nelems)
{
  assert(nelems);
  line_t line;
  seek_prefix(f, line, sizeof(line), "<Cells");
  line_t name;
  while (1) {
    seek_prefix(f, line, sizeof(line), "<DataArray");
    read_array_name(line, name);
    if (!strcmp(name, "types"))
      break;
  }
  enum tag_type type;
  unsigned ncomps;
  void* data;
  read_array(f, line, &type, name, nelems, &ncomps, &data);
  assert(type == TAG_U8);
  assert(!strcmp(name, "types"));
  assert(ncomps == 1);
  unsigned char* types = (unsigned char*) data;
  unsigned dim;
  for (dim = 0; dim < 4; ++dim)
    if (types[0] == simplex_types[dim])
      break;
  assert(dim < 4);
  for (unsigned i = 1; i < nelems; ++i)
    assert(types[i] == simplex_types[dim]);
  loop_host_free(types);
  return dim;
}

static unsigned read_tags(FILE* f, char const* prefix, struct tags* ts,
    unsigned n)
{
  line_t line;
  seek_prefix(f, line, sizeof(line), prefix);
  unsigned nt = 0;
  while(read_tag(f, ts, n))
    ++nt;
  return nt;
}

static void read_points(FILE* f, struct tags* ts, unsigned n)
{
  unsigned nt = read_tags(f, "<Points", ts, n);
  assert(nt == 1);
}

static void read_verts(FILE* f, struct mesh* m)
{
  read_points(f, mesh_tags(m, 0), mesh_count(m, 0));
}

static char const* get_dim_name(struct mesh* m, unsigned ent_dim)
{
  if (mesh_dim(m) == ent_dim)
    return "Cell";
  else
    switch (ent_dim) {
      case 0:
        return "Point";
      case 1:
        return "Edge";
      case 2:
        return "Face";
      default:
        return 0;
    }
}

static void read_ents(FILE* f, struct mesh* m, unsigned nents,
    unsigned ent_dim)
{
  line_t tag;
  sprintf(tag, "<%s", get_dim_name(m, ent_dim));
  line_t line;
  seek_prefix(f, line, sizeof(line), tag);
  seek_prefix(f, line, sizeof(line), "<DataArray");
  unsigned verts_per_elem = the_down_degrees[ent_dim][0];
  enum tag_type type;
  line_t name;
  unsigned ncomps;
  void* data;
  read_array(f, line, &type, name, nents * verts_per_elem, &ncomps, &data);
  assert(type == TAG_U32);
  assert(!strcmp("connectivity", name));
  assert(ncomps == 1);
  mesh_set_ents(m, ent_dim, nents, (unsigned*) data);
}

static struct mesh* read_vtk_mesh(FILE* f)
{
  unsigned nverts, nelems;
  unsigned do_edges, do_faces;
  unsigned nedges, nfaces;
  read_ent_counts(f, &nverts, &nelems,
      &do_edges, &do_faces, &nedges, &nfaces);
  assert(nelems);
  unsigned dim = read_dimension(f, nelems);
  struct mesh* m = new_mesh(dim);
  mesh_set_ents(m, 0, nverts, 0);
  rewind(f);
  read_verts(f, m);
  if (do_edges)
    read_ents(f, m, nedges, 1);
  if (do_faces)
    read_ents(f, m, nfaces, 2);
  read_ents(f, m, nelems, dim);
  return m;
}

static void read_vtk_fields(FILE* f, struct mesh* m)
{
  unsigned dim = mesh_dim(m);
  read_tags(f, "<PointData", mesh_tags(m, 0), mesh_count(m, 0));
  if ((dim > 1) && mesh_has_dim(m, 1))
    read_tags(f, "<EdgeData", mesh_tags(m, 1), mesh_count(m, 1));
  if ((dim > 2) && mesh_has_dim(m, 2))
    read_tags(f, "<FaceData", mesh_tags(m, 2), mesh_count(m, 2));
  read_tags(f, "<CellData", mesh_tags(m, dim), mesh_count(m, dim));
}

struct mesh* read_vtu(char const* filename)
{
  FILE* f = fopen(filename, "r");
  assert(f);
  struct mesh* m = read_vtk_mesh(f);
  read_vtk_fields(f, m);
  fclose(f);
  return m;
}

void write_vtu_cloud_opts(struct cloud* c, char const* filename,
    enum vtk_format fmt)
{
  unsigned npts = cloud_count(c);
  FILE* file = fopen(filename, "w");
  write_unstructured_header(file);
  fprintf(file, "<UnstructuredGrid>\n");
  fprintf(file, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"1\">\n", npts);
  fprintf(file, "<Points>\n");
  struct const_tag* coord_tag = cloud_find_tag(c, "coordinates");
  write_tag(file, npts, coord_tag, fmt);
  fprintf(file, "</Points>\n");
  fprintf(file, "<Cells>\n");
  unsigned* conn = LOOP_HOST_MALLOC(unsigned, npts);
  for (unsigned i = 0; i < npts; ++i)
    conn[i] = i;
  write_array(file, TAG_U32, "connectivity", npts, 1, conn, fmt);
  loop_host_free(conn);
  unsigned off[1];
  off[0] = npts;
  write_array(file, TAG_U32, "offsets", 1, 1, off, fmt);
  unsigned char type[1] = {VTK_POLY_VERTEX};
  write_array(file, TAG_U8, "types", 1, 1, type, fmt);
  fprintf(file, "</Cells>\n");
  fprintf(file, "<PointData>\n");
  for (unsigned i = 0; i < cloud_count_tags(c); ++i) {
    struct const_tag* tag = cloud_get_tag(c, i);
    if (tag != coord_tag)
      write_tag(file, npts, tag, fmt);
  }
  fprintf(file, "</PointData>\n");
  fprintf(file, "<CellData>\n");
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n");
  fprintf(file, "</UnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");
  fclose(file);
}

void write_vtu_cloud(struct cloud* c, char const* filename)
{
  write_vtu_cloud_opts(c, filename, VTK_BINARY);
}

struct cloud* read_vtu_cloud(char const* filename)
{
  FILE* f = fopen(filename, "r");
  assert(f);
  unsigned npts;
  read_nverts(f, &npts);
  assert(npts);
  struct cloud* c = new_cloud(npts);
  read_points(f, cloud_tags(c), npts);
  read_tags(f, "<PointData", cloud_tags(c), npts);
  fclose(f);
  return c;
}

static void write_pieces(FILE* file, char const* pathname, unsigned npieces)
{
  line_t a;
  char* filename;
  split_pathname(pathname, a, sizeof(a), &filename, 0);
  for (unsigned i = 0; i < npieces; ++i) {
    line_t b;
    enum_pathname(filename, npieces, i, "vtu", b, sizeof(b));
    fprintf(file, "<Piece Source=\"%s\"/>\n", b);
  }
}

void write_pvtu(struct mesh* m, char const* filename,
    unsigned npieces)
{
  FILE* file = fopen(filename, "w");
  fprintf(file, "<VTKFile type=\"PUnstructuredGrid\">\n");
  fprintf(file, "<PUnstructuredGrid GhostLevel=\"%u\">\n", mesh_ghost_layers(m));
  struct const_tag* coord_tag = mesh_find_tag(m, 0, "coordinates");
  fprintf(file, "<PPointData>\n");
  for (unsigned i = 0; i < mesh_count_tags(m, 0); ++i) {
    struct const_tag* t = mesh_get_tag(m, 0, i);
    if (t != coord_tag) {
      fprintf(file, "<PDataArray ");
      describe_tag(file, t);
      fprintf(file, "/>\n");
    }
  }
  fprintf(file, "</PPointData>\n");
  fprintf(file, "<PCellData>\n");
  for (unsigned i = 0; i < mesh_count_tags(m, mesh_dim(m)); ++i) {
    fprintf(file, "<PDataArray ");
    describe_tag(file, mesh_get_tag(m, mesh_dim(m), i));
    fprintf(file, "/>\n");
  }
  fprintf(file, "</PCellData>\n");
  fprintf(file, "<PPoints>\n");
  fprintf(file, "<PDataArray ");
  describe_tag(file, coord_tag);
  fprintf(file, "/>\n");
  fprintf(file, "</PPoints>\n");
  write_pieces(file, filename, npieces);
  fprintf(file, "</PUnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");
  fclose(file);
}

void write_pvtu_cloud(struct cloud* c, char const* filename,
    unsigned npieces)
{
  FILE* file = fopen(filename, "w");
  fprintf(file, "<VTKFile type=\"PUnstructuredGrid\">\n");
  fprintf(file, "<PUnstructuredGrid>\n");
  struct const_tag* coord_tag = cloud_find_tag(c, "coordinates");
  fprintf(file, "<PPointData>\n");
  for (unsigned i = 0; i < cloud_count_tags(c); ++i) {
    struct const_tag* t = cloud_get_tag(c, i);
    if (t != coord_tag) {
      fprintf(file, "<PDataArray ");
      describe_tag(file, t);
      fprintf(file, "/>\n");
    }
  }
  fprintf(file, "</PPointData>\n");
  fprintf(file, "<PPoints>\n");
  fprintf(file, "<PDataArray ");
  describe_tag(file, coord_tag);
  fprintf(file, "/>\n");
  fprintf(file, "</PPoints>\n");
  write_pieces(file, filename, npieces);
  fprintf(file, "</PUnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");
  fclose(file);
}

struct mesh* read_parallel_vtu(char const* inpath)
{
  char* suffix;
  line_t prefix;
  split_pathname(inpath, prefix, sizeof(prefix), 0, &suffix);
  line_t piecepath;
  enum_pathname(prefix, comm_size(), comm_rank(), "vtu",
      piecepath, sizeof(piecepath));
  struct mesh* m = read_vtu(piecepath);
  if (mesh_find_tag(m, mesh_dim(m), "piece"))
    mesh_free_tag(m, mesh_dim(m), "piece");
  return m;
}

void write_parallel_vtu(struct mesh* m, char const* outpath)
{
  char* suffix;
  line_t prefix;
  split_pathname(outpath, prefix, sizeof(prefix), 0, &suffix);
  line_t piecepath;
  enum_pathname(prefix, comm_size(), comm_rank(), "vtu",
      piecepath, sizeof(piecepath));
  unsigned* piece = uints_copy(
      mesh_ask_own_ranks(m, mesh_dim(m)),
      mesh_count(m, mesh_dim(m)));
  mesh_add_tag(m, mesh_dim(m), TAG_U32, "piece", 1, piece);
  write_vtu(m, piecepath);
  if (!comm_rank() && !strcmp(suffix, "pvtu"))
    write_pvtu(m, outpath, comm_size());
  mesh_free_tag(m, mesh_dim(m), "piece");
}

struct cloud* read_parallel_vtu_cloud(char const* inpath)
{
  char* suffix;
  line_t prefix;
  split_pathname(inpath, prefix, sizeof(prefix), 0, &suffix);
  line_t piecepath;
  enum_pathname(prefix, comm_size(), comm_rank(), "vtu",
      piecepath, sizeof(piecepath));
  struct cloud* c = read_vtu_cloud(piecepath);
  if (cloud_find_tag(c, "piece"))
    cloud_free_tag(c, "piece");
  return c;
}

void write_parallel_vtu_cloud(struct cloud* c, char const* outpath)
{
  char* suffix;
  line_t prefix;
  split_pathname(outpath, prefix, sizeof(prefix), 0, &suffix);
  line_t piecepath;
  enum_pathname(prefix, comm_size(), comm_rank(), "vtu",
      piecepath, sizeof(piecepath));
  unsigned* piece = uints_filled(cloud_count(c), comm_rank());
  cloud_add_tag(c, TAG_U32, "piece", 1, piece);
  write_vtu_cloud(c, piecepath);
  if (!comm_rank() && !strcmp(suffix, "pvtu"))
    write_pvtu_cloud(c, outpath, comm_size());
  cloud_free_tag(c, "piece");
}
