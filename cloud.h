#ifndef CLOUD_H
#define CLOUD_H

#include "tag.h"

struct cloud* new_cloud(unsigned count);
void free_cloud(struct cloud* c);

unsigned cloud_count(struct cloud* c);

struct const_tag* cloud_add_tag(struct cloud* c, enum tag_type type,
    char const* name, unsigned ncomps, void* data);
void cloud_free_tag(struct cloud* c, char const* name);
struct const_tag* cloud_find_tag(struct cloud* c, char const* name);
unsigned cloud_count_tags(struct cloud* c);
struct const_tag* cloud_get_tag(struct cloud* c, unsigned i);
void cloud_rename_tag(struct cloud* c, char const* oldname,
    char const* newname);

#endif
