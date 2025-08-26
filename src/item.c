#include "item.h"
#include <math.h>

void fill_list(Item *buff, const Instance *restrict const inst) {
    for (size_t i = 0; i < inst->n_items; i++) {
        buff->a_val = inst->a_ptr[i];
        buff->b_val = inst->b_ptr[i];
        buff->p_val = inst->p_ptr[i];
        buff->id = i;
        buff->selected = false;
        buff++;
    }
}

typedef int (*Sorter)(const void *, const void *);

int get_order(double a, double b) {
    if (a < b)
        return 1;
    if (b < a)
        return -1;
    return 0;
}

int a_weight_sort(const void *a, const void *b) {
    const Item *a_ptr = a;
    const Item *b_ptr = b;

    double a_cost = a_ptr->p_val / a_ptr->a_val;
    double b_cost = b_ptr->p_val / b_ptr->a_val;
    //double a_cost = a_ptr->p_val;
    //double b_cost = b_ptr->p_val;
    return get_order(a_cost, b_cost);
}

int full_weight_sort(const void *a, const void *b) {
    const Item *a_ptr = a;
    const Item *b_ptr = b;

    double a_cost = a_ptr->p_val / (a_ptr->a_val + sqrt(a_ptr->b_val));
    double b_cost = b_ptr->p_val / (b_ptr->a_val + sqrt(b_ptr->b_val));
    //double a_cost = a_ptr->p_val;
    //double b_cost = b_ptr->p_val;
    return get_order(a_cost, b_cost);
}

Sorter get_sorter(const ItemSort sort) {

    Sorter sorter = NULL;
    switch (sort) {
    case AWeight:
        sorter = a_weight_sort;
        break;
    case FullWeight:
        sorter = full_weight_sort;
        break;
    case None:
        sorter = NULL;
        break;
    }

    return sorter;
}

void item_list_sort(ItemList *item_list, const ItemSort sort) {
    Sorter sorter = get_sorter(sort);
    qsort(item_list->item_ptr, item_list->n_items, sizeof(Item), sorter);
}

int item_list_alloc(ItemList *item_list,
                    const Instance *restrict const instance,
                    const ItemSort sort) {

    Item *buff = malloc(sizeof(Item) * instance->n_items);
    if (buff == NULL)
        return FAILURE;

    item_list->item_ptr = buff;
    fill_list(buff, instance);
    item_list->n_items = instance->n_items;
    if (sort != None) 
        item_list_sort(item_list, sort);
   


    return SUCCESS;
}

void item_list_free(ItemList *item_list) {
    free(item_list->item_ptr);
    item_list->item_ptr = NULL;
}

void invert_item_list(ItemList *item_list) {

    size_t begin = 0;
    size_t end = item_list->n_items - 1;
    while (begin < end) {
        Item tmp = item_list->item_ptr[begin];
        item_list->item_ptr[begin] = item_list->item_ptr[end];
        item_list->item_ptr[end] = tmp;
        begin++;
        end--;
    }
}
