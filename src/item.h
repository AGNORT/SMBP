#ifndef ITEM_H
#define ITEM_H
#ifndef RESTRICT
#  if defined(__cplusplus)
#    define RESTRICT
#  else
#    define RESTRICT restrict
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "instance.h"
#include <stdbool.h>

typedef struct {
    size_t id;
    double a_val;
    double b_val;
    double p_val;
    bool selected;
} Item;

typedef struct {
    Item *item_ptr;
    size_t n_items;
} ItemList;

typedef enum { None, AWeight, FullWeight } ItemSort;

int item_list_alloc(ItemList *item_list,
                    const Instance *RESTRICT const instance,
                    const ItemSort sort);

void item_list_sort(ItemList *item_list, const ItemSort sort);
void invert_item_list(ItemList *item_list);

void item_list_free(ItemList *item_list);

#ifdef __cplusplus
}
#endif
#endif
