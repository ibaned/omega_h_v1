#ifndef CLOUD_H
#define CLOUD_H

struct cloud* new_cloud(unsigned count);
void free_cloud(struct cloud* c);

unsigned cloud_count(struct cloud* c);

struct const_field* cloud_add_field(struct cloud* c, char const* name,
    unsigned ncomps, double* data);
void cloud_free_field(struct cloud* c, char const* name);
struct const_field* cloud_find_field(struct cloud* c, char const* name);
unsigned cloud_count_fields(struct cloud* c);
struct const_field* cloud_get_field(struct cloud* c, unsigned i);

struct const_label* cloud_add_label(struct cloud* c, char const* name,
    unsigned* data);
void cloud_free_label(struct cloud* c, char const* name);
struct const_label* cloud_find_label(struct cloud* c, char const* name);
unsigned cloud_count_labels(struct cloud* c);
struct const_label* cloud_get_label(struct cloud* c, unsigned i);

#endif
