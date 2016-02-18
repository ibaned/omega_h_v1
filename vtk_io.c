#include "vtk_io.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arrays.h"
#include "base64.h"
#include "comm.h"
#include "compress.h"
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
  }
  LOOP_NORETURN("");
}

static char const* format_name(enum vtk_format fmt)
{
  switch (fmt) {
    case VTK_ASCII: return "ascii";
    case VTK_BINARY: return "binary";
  }
  LOOP_NORETURN("");
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
  LOOP_NORETURN(TAG_U32);
}

static enum vtk_format read_array_format(char const* header)
{
  line_t val;
  read_attrib(header, "format", val);
  for (unsigned fmt = 0; fmt < VTK_FORMATS; ++fmt)
    if (!strcmp(format_name((enum vtk_format) fmt), val))
      return (enum vtk_format) fmt;
  LOOP_NORETURN(VTK_ASCII);
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

static void write_binary_uints(FILE* file, unsigned const* x, unsigned n)
{
  char* s = base64_encode(x, n * sizeof(unsigned));
  fputs(s, file);
  loop_host_free(s);
}

static void write_binary_array(FILE* file, enum tag_type t, unsigned nents,
    unsigned ncomps, void const* data)
{
  unsigned tsize = tag_size(t);
  unsigned size = tsize * ncomps * nents;
  unsigned long comp_size;
  void* host_data = uchars_to_host((unsigned char const*) data, size);
  void* comp = my_compress(host_data, size, &comp_size);
  loop_host_free(host_data);
  if (can_compress) {
    unsigned comp_header[4] = {1, size, size, (unsigned) comp_size};
    write_binary_uints(file, comp_header, 4);
  } else
    write_binary_uints(file, &size, 1);
  char* s = base64_encode(comp, comp_size);
  loop_host_free(comp);
  fputs(s, file);
  loop_host_free(s);
  fputc('\n', file);
}

static void read_binary_uints(char const** p, unsigned* x, unsigned n)
{
  unsigned* dec = (unsigned*) base64_decode(p, n * sizeof(unsigned));
  for (unsigned i = 0; i < n; ++i)
    x[i] = dec[i];
  loop_host_free(dec);
}

static void* read_binary_array(FILE* file, enum endian end, unsigned do_com,
    enum tag_type t, unsigned nents, unsigned ncomps)
{
  unsigned tsize = tag_size(t);
  unsigned long enc_nchars;
  char* enc = base64_fread(file, &enc_nchars);
  char const* p = enc;
  unsigned long decomp_size = nents * ncomps * tsize;
  unsigned long comp_size;
  if (do_com) {
    unsigned comp_header[4];
    read_binary_uints(&p, comp_header, 4);
    if (end != endianness())
      swap_one(&comp_header[3], sizeof(unsigned));
    comp_size = comp_header[3];
  } else {
    unsigned ignore;
    read_binary_uints(&p, &ignore, 1);
    comp_size = decomp_size;
  }
  void* decod = base64_decode(&p, comp_size);
  loop_host_free(enc);
  void* decomp;
  if (do_com) {
    decomp = my_decompress(decod, comp_size, decomp_size);
    loop_host_free(decod);
  } else
    decomp = decod;
  /* this function also copies to device space */
  void* swapped = generic_swap_if_needed(end, nents * ncomps, tsize, decomp);
  loop_host_free(decomp);
  return swapped;
}

static void write_ascii_array(FILE* file, enum tag_type t, unsigned nents,
    unsigned ncomps, void const* data)
{
  switch (t) {
    case TAG_U8: {
      unsigned char* p = uchars_to_host((unsigned char const*) data, nents * ncomps);
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < ncomps; ++j)
          fprintf(file, " %hhu", p[i * ncomps + j]);
        fprintf(file, "\n");
      }
      loop_host_free(p);
      break;
    }
    case TAG_U32: {
      unsigned* p = uints_to_host((unsigned int const*) data, nents * ncomps);
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < ncomps; ++j)
          fprintf(file, " %u", p[i * ncomps + j]);
        fprintf(file, "\n");
      }
      loop_host_free(p);
      break;
    }
    case TAG_U64: {
      unsigned long* p = ulongs_to_host((unsigned long const*) data, nents * ncomps);
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < ncomps; ++j)
          fprintf(file, " %lu", p[i * ncomps + j]);
        fprintf(file, "\n");
      }
      loop_host_free(p);
      break;
    }
    case TAG_F64: {
      double* p = doubles_to_host((double const*) data, nents * ncomps);
      for (unsigned i = 0; i < nents; ++i) {
        for (unsigned j = 0; j < ncomps; ++j)
          fprintf(file, " %.15e", p[i * ncomps + j]);
        fprintf(file, "\n");
      }
      loop_host_free(p);
      break;
    }
  }
}

static void* read_ascii_array(FILE* file, enum tag_type type, unsigned nents,
    unsigned ncomps)
{
  unsigned n = nents * ncomps;
  void* out = 0;
  switch (type) {
    case TAG_U8: {
      unsigned char* in = LOOP_HOST_MALLOC(unsigned char, n);
      for (unsigned i = 0; i < n; ++i)
        safe_scanf(file, 1, "%hhu", &in[i]);
      out = uchars_to_device(in, n);
      loop_host_free(in);
      break;
    }
    case TAG_U32: {
      unsigned* in = LOOP_HOST_MALLOC(unsigned, n);
      for (unsigned i = 0; i < n; ++i)
        safe_scanf(file, 1, "%u", &in[i]);
      out = uints_to_device(in, n);
      loop_host_free(in);
      break;
    }
    case TAG_U64: {
      unsigned long* in = LOOP_HOST_MALLOC(unsigned long, n);
      for (unsigned i = 0; i < n; ++i)
        safe_scanf(file, 1, "%lu", &in[i]);
      out = ulongs_to_device(in, n);
      loop_host_free(in);
      break;
    }
    case TAG_F64: {
      double* in = LOOP_HOST_MALLOC(double, n);
      for (unsigned i = 0; i < n; ++i)
        safe_scanf(file, 1, "%lf", &in[i]);
      out = doubles_to_device(in, n);
      loop_host_free(in);
      break;
    }
  }
  return out;
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

static void read_array(FILE* f, char const* line,
    enum endian end, unsigned do_com,
    enum tag_type* type, char* name,
    unsigned nents,
    unsigned* ncomps, void** data)
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
      *data = read_binary_array(f, end, do_com, *type, nents, *ncomps);
      break;
  }
  line_t tmpline;
  seek_prefix(f, tmpline, sizeof(tmpline), "</DataArray");
}

static unsigned read_tag(FILE* f, struct tags* ts, unsigned n,
    enum endian end, unsigned do_com)
{
  line_t line;
  if (!seek_prefix_next(f, line, sizeof(line), "<DataArray"))
    return 0;
  enum tag_type type;
  line_t name;
  unsigned ncomps;
  void* data;
  read_array(f, line, end, do_com, &type, name, n, &ncomps, &data);
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
  unsigned nverts_per_elem = the_down_degrees[elem_dim][0];
  unsigned* off = uints_linear(nelems + 1, nverts_per_elem);
  write_array(file, TAG_U32, "offsets", nelems, 1, off + 1, fmt);
  loop_free(off);
  unsigned char type = (unsigned char) simplex_types[elem_dim];
  unsigned char* types = uchars_filled(nelems, type);
  write_array(file, TAG_U8, "types", nelems, 1, types, fmt);
  loop_free(types);
}

static void write_mesh_tags(FILE* file, struct mesh* m, unsigned dim,
    enum vtk_format fmt, struct const_tag* except)
{
  for (unsigned i = 0; i < mesh_count_tags(m, dim); ++i) {
    struct const_tag* tag = mesh_get_tag(m, dim, i);
    if (tag != except)
      write_tag(file, mesh_count(m, dim), tag, fmt);
  }
}

static void write_unstructured_header(FILE* file, enum vtk_format fmt)
{
  fprintf(file, "<VTKFile type=\"UnstructuredGrid\"");
  if (fmt == VTK_BINARY) {
    if (endianness() == MY_LITTLE_ENDIAN)
      fprintf(file, " byte_order=\"LittleEndian\"");
    else
      fprintf(file, " byte_order=\"BigEndian\"");
    fprintf(file, " header_type=\"UInt32\"");
    if (can_compress)
      fprintf(file, " compressor=\"vtkZLibDataCompressor\"");
  }
  fprintf(file, ">\n");
}

static void read_unstructured_header(FILE* file,
    enum endian* end, unsigned* do_com)
{
  line_t line;
  seek_prefix(file, line, sizeof(line), "<VTKFile");
  line_t byte_order;
  if (try_read_attrib(line, "byte_order", byte_order)) {
    if (!strcmp(byte_order, "LittleEndian"))
      *end = MY_LITTLE_ENDIAN;
    else
      *end = MY_BIG_ENDIAN;
  } else
    *end = endianness();
  line_t compressor;
  *do_com = try_read_attrib(line, "compressor", compressor);
}

static void write_piece_header(FILE* file, struct mesh* m)
{
  unsigned elem_dim = mesh_dim(m);
  fprintf(file, "<Piece NumberOfPoints=\"%u\" NumberOfCells=\"%u\"",
      mesh_count(m, 0), mesh_count(m, elem_dim));
  if ((elem_dim > 1) && mesh_has_dim(m, 1))
    fprintf(file, " NumberOfEdges=\"%u\"", mesh_count(m, 1));
  if ((elem_dim > 2) && mesh_has_dim(m, 2))
    fprintf(file, " NumberOfFaces=\"%u\"", mesh_count(m, 2));
  fprintf(file, " Rep=\"%s\"",
      mesh_get_rep(m) == MESH_FULL ? "Full" : "Reduced");
  fprintf(file, ">\n");
}

static void read_piece_header(FILE* f, unsigned* nverts, unsigned* nelems,
    unsigned* do_edges, unsigned* do_faces,
    unsigned* nedges, unsigned* nfaces,
    enum mesh_rep* rep)
{
  line_t line;
  seek_prefix(f, line, sizeof(line), "<Piece");
  *nverts = read_int_attrib(line, "NumberOfPoints");
  *nelems = read_int_attrib(line, "NumberOfCells");
  *do_edges = try_read_int_attrib(line, "NumberOfEdges", nedges);
  *do_faces = try_read_int_attrib(line, "NumberOfFaces", nfaces);
  line_t rep_text;
  if (try_read_attrib(line, "Rep", rep_text) &&
      !strcmp(rep_text, "Full"))
    *rep = MESH_FULL;
  else
    *rep = MESH_REDUCED;
}

static void write_vtu_opts(struct mesh* m, char const* filename, enum vtk_format fmt)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nverts = mesh_count(m, 0);
  unsigned do_edges = ((elem_dim > 1) && mesh_has_dim(m, 1));
  unsigned do_faces = ((elem_dim > 2) && mesh_has_dim(m, 2));
  FILE* file = safe_fopen(filename, "w");
  write_unstructured_header(file, fmt);
  fprintf(file, "<UnstructuredGrid>\n");
  write_piece_header(file, m);
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
  write_mesh_tags(file, m, 0, fmt, coord_tag);
  fprintf(file, "</PointData>\n");
  if (do_edges) {
    fprintf(file, "<EdgeData>\n");
    write_mesh_tags(file, m, 1, fmt, 0);
    fprintf(file, "</EdgeData>\n");
  }
  if (do_faces) {
    fprintf(file, "<FaceData>\n");
    write_mesh_tags(file, m, 2, fmt, 0);
    fprintf(file, "</FaceData>\n");
  }
  fprintf(file, "<CellData>\n");
  write_mesh_tags(file, m, elem_dim, fmt, 0);
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n");
  fprintf(file, "</UnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");
  fclose(file);
}

static void write_vtu(struct mesh* m, char const* filename)
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

LOOP_KERNEL(assert_one_type,
    unsigned char* types,
    unsigned char type)
  assert(types[i] == type);
}

static void seek_types(FILE* f, line_t line, line_t name)
{
  seek_prefix(f, line, sizeof(line_t), "<Cells");
  while (1) {
    seek_prefix(f, line, sizeof(line_t), "<DataArray");
    read_array_name(line, name);
    if (!strcmp(name, "types"))
      return;
  }
}

static unsigned read_dimension(FILE* f, unsigned nelems, enum endian end,
    unsigned do_com)
{
  assert(nelems);
  line_t line;
  line_t name;
  seek_types(f, line, name);
  enum tag_type type;
  unsigned ncomps;
  void* data;
  read_array(f, line, end, do_com, &type, name, nelems, &ncomps, &data);
  assert(type == TAG_U8);
  assert(!strcmp(name, "types"));
  assert(ncomps == 1);
  unsigned char* types = (unsigned char*) data;
  unsigned dim;
  unsigned char first_type = uchars_at(types, 0);
  for (dim = 0; dim < 4; ++dim)
    if (first_type == simplex_types[dim])
      break;
  assert(dim < 4);
  LOOP_EXEC(assert_one_type, nelems, types, simplex_types[dim]);
  loop_free(types);
  return dim;
}

static unsigned read_tags(FILE* f, char const* prefix, struct tags* ts,
    unsigned n, enum endian end, unsigned do_com)
{
  line_t line;
  seek_prefix(f, line, sizeof(line), prefix);
  unsigned nt = 0;
  while(read_tag(f, ts, n, end, do_com))
    ++nt;
  return nt;
}

static void read_points(FILE* f, struct tags* ts, unsigned n, enum endian end,
    unsigned do_com)
{
  unsigned nt = read_tags(f, "<Points", ts, n, end, do_com);
  assert(nt == 1);
}

static void read_verts(FILE* f, struct mesh* m, enum endian end,
    unsigned do_com)
{
  read_points(f, mesh_tags(m, 0), mesh_count(m, 0), end, do_com);
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
    unsigned ent_dim, enum endian end, unsigned do_com)
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
  read_array(f, line, end, do_com,
      &type, name, nents * verts_per_elem, &ncomps, &data);
  assert(type == TAG_U32);
  assert(!strcmp("connectivity", name));
  assert(ncomps == 1);
  mesh_set_ents(m, ent_dim, nents, (unsigned*) data);
}

static struct mesh* read_vtk_mesh(FILE* f, enum endian end,
    unsigned do_com, unsigned is_parallel)
{
  unsigned nverts, nelems;
  unsigned do_edges, do_faces;
  unsigned nedges, nfaces;
  enum mesh_rep rep;
  read_piece_header(f, &nverts, &nelems,
      &do_edges, &do_faces, &nedges, &nfaces,
      &rep);
  assert(nelems);
  unsigned dim = read_dimension(f, nelems, end, do_com);
  struct mesh* m = new_mesh(dim, rep, is_parallel);
  mesh_set_ents(m, 0, nverts, 0);
  rewind(f);
  read_verts(f, m, end, do_com);
  if (do_edges)
    read_ents(f, m, nedges, 1, end, do_com);
  if (do_faces)
    read_ents(f, m, nfaces, 2, end, do_com);
  read_ents(f, m, nelems, dim, end, do_com);
  return m;
}

static void read_mesh_tags(FILE* f, struct mesh* m, unsigned dim,
    char const* prefix, enum endian end, unsigned do_com, unsigned is_parallel)
{
  read_tags(f, prefix, mesh_tags(m, dim), mesh_count(m, dim), end, do_com);
  if (is_parallel)
    mesh_parallel_from_tags(m, dim);
}

static void read_vtk_fields(FILE* f, struct mesh* m, enum endian end,
    unsigned do_com, unsigned is_parallel)
{
  unsigned dim = mesh_dim(m);
  read_mesh_tags(f, m, 0, "<PointData", end, do_com, is_parallel);
  if ((dim > 1) && mesh_has_dim(m, 1))
    read_mesh_tags(f, m, 1, "<EdgeData", end, do_com, is_parallel);
  if ((dim > 2) && mesh_has_dim(m, 2))
    read_mesh_tags(f, m, 2, "<FaceData", end, do_com, is_parallel);
  read_mesh_tags(f, m, dim, "<CellData", end, do_com, is_parallel);
}

static struct mesh* read_vtu_opts(char const* filename, unsigned is_parallel)
{
  FILE* file = safe_fopen(filename, "r");
  enum endian end;
  unsigned do_com;
  read_unstructured_header(file, &end, &do_com);
  if (do_com)
    assert(can_compress);
  struct mesh* m = read_vtk_mesh(file, end, do_com, is_parallel);
  read_vtk_fields(file, m, end, do_com, is_parallel);
  fclose(file);
  return m;
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

static void write_pvtu(struct mesh* m, char const* filename,
    unsigned npieces)
{
  FILE* file = safe_fopen(filename, "w");
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

static void get_ghost_level(struct mesh* m, char const* pvtu_filename)
{
  FILE* file = safe_fopen(pvtu_filename, "r");
  line_t line;
  seek_prefix(file, line, sizeof(line), "<PUnstructuredGrid");
  unsigned nlayers;
  if (try_read_int_attrib(line, "GhostLevel", &nlayers))
    mesh_set_ghost_layers(m, nlayers);
  fclose(file);
}

struct mesh* read_mesh_vtk(char const* inpath)
{
  char* suffix;
  line_t prefix;
  split_pathname(inpath, prefix, sizeof(prefix), 0, &suffix);
  unsigned is_parallel = !(strcmp(suffix, "pvtu"));
  if (!is_parallel && comm_size() > 1)
    /* fine, but I'm not going to like it. */
    return read_and_partition_serial_mesh(inpath);
  line_t piecepath;
  enum_pathname(prefix, comm_size(), comm_rank(), "vtu",
      piecepath, sizeof(piecepath));
  struct mesh* m = read_vtu_opts(piecepath, is_parallel);
  if (is_parallel)
    get_ghost_level(m, inpath);
  return m;
}

void write_mesh_vtk_opts(struct mesh* m, char const* outpath,
    enum vtk_format fmt)
{
  char* suffix;
  line_t prefix;
  split_pathname(outpath, prefix, sizeof(prefix), 0, &suffix);
  line_t piecepath;
  enum_pathname(prefix, comm_size(), comm_rank(), "vtu",
      piecepath, sizeof(piecepath));
  if (mesh_is_parallel(m))
    for (unsigned dim = 0; dim <= mesh_dim(m); ++dim)
      mesh_parallel_to_tags(m, dim);
  write_vtu_opts(m, piecepath, fmt);
  if (!comm_rank() && !strcmp(suffix, "pvtu"))
    write_pvtu(m, outpath, comm_size());
  if (mesh_is_parallel(m))
    for (unsigned dim = 0; dim <= mesh_dim(m); ++dim)
      mesh_parallel_untag(m, dim);
}

void write_mesh_vtk(struct mesh* m, char const* outpath)
{
  write_mesh_vtk_opts(m, outpath, VTK_BINARY);
}
