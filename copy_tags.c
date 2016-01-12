#include "copy_tags.h"

#include "arrays.h"
#include "tag.h"

void copy_tags(struct tags* a, struct tags* b, unsigned n)
{
  for (unsigned i = 0; i < count_tags(a); ++i) {
    struct const_tag* t = get_tag(a, i);
    void* data = 0;
    switch (t->type) {
      case TAG_U8:  data = uchars_copy(t->d.u8, n * t->ncomps);
                    break;
      case TAG_U32: data = uints_copy(t->d.u32, n * t->ncomps);
                    break;
      case TAG_U64: data = ulongs_copy(t->d.u64, n * t->ncomps);
                    break;
      case TAG_F64: data = doubles_copy(t->d.f64, n * t->ncomps);
                    break;
    };
    add_tag(b, t->type, t->name, t->ncomps, data);
  }
}
