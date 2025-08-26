#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ts_config.h"
#include "instance.h"

#define INF 1e100

int int_in_range(int lb, int ub) {
    int rnd = rand();
    int bnd = ub - lb;
    if (bnd == 0)
        return lb;
    return (rnd % bnd) + lb;
}

typedef struct {
    bool *items;
    double obj_val;
    double a_weight;
    double b_weight;
    size_t n_items;
} Solution;

void solution_free(Solution *restrict const sol) { free(sol->items); }

int solution_init(Solution *restrict sol, size_t n_items) {

    bool *buff = calloc(n_items, sizeof(bool));
    if (buff == NULL)
        return 1;

    sol->items = buff;
    sol->obj_val = 0;
    sol->n_items = n_items;
    return 0;
}

int solution_init_from(const Solution *restrict const src,
                       Solution *restrict dst) {

    dst->obj_val = src->obj_val;
    dst->n_items = src->n_items;
    bool *buff = calloc(src->n_items, sizeof(bool));
    if (buff == NULL)
        return 1;

    for (size_t i = 0; i < src->n_items; i++) {
        buff[i] = src->items[i];
    }
    dst->items = buff;

    return 0;
}

void initial_solution(Solution *restrict sol,
                      const Instance *const restrict inst) {

    double tot_a = 0;
    double tot_b = 0;

    for (size_t i = 0; i < inst->n_items; i++) {

        double curr_a = tot_a + inst->a_ptr[i];
        double curr_b = tot_b + inst->b_ptr[i];
        double curr_w = curr_a + inst->rho * sqrt(curr_b);
        if (curr_w <= inst->capacity) {
            sol->items[i] = true;
            sol->obj_val += inst->p_ptr[i];
            tot_a = curr_a;
            tot_b = curr_b;
        } else {
            break;
        }
    }

    sol->a_weight = tot_a;
    sol->b_weight = tot_b;
}

typedef struct {
    int *block_move;
    size_t n_items;
    int block_size;
} TabuTenure;

void tt_free(TabuTenure *const restrict tt) { free(tt->block_move); }

void tt_reset(TabuTenure *tt) {
    for (size_t i = 0; i < tt->n_items; i++)
        tt->block_move[i] = 0;
}

int max(int a, int b) { return (a > b) ? a : b; }

int tabu_tenure_init(TabuTenure *restrict tt, size_t n_items) {
    tt->n_items = n_items;
    int *buff = calloc(n_items, sizeof(int));
    if (buff == NULL)
        return 1;
    tt->block_move = buff;
    tt->block_size = max(n_items / BLOCK_SCALE, MIN_BLOCK_SZ);

    return 0;
}

void tt_update(TabuTenure *tt) {
    int *tmp = tt->block_move;
    size_t count = tt->n_items;
    bool need_reset = true;
    while (count--) {
        if (*tmp)
            (*tmp)--;
        else
            need_reset = false;
        tmp++;
    }

    if (need_reset) {
    }
}

typedef struct {
    Solution *best;
    Solution *curr;
} WorkSols;

void work_sols_swap(WorkSols *work_sol) {
    Solution *tmp = work_sol->best;
    work_sol->best = work_sol->curr;
    work_sol->curr = tmp;
}

bool tt_is_free(const TabuTenure *restrict const tt, int item) {
    return tt->block_move[item] == 0;
}

double item_weight(const Instance *const restrict inst, int idx) {
    return inst->a_ptr[idx] + inst->rho * inst->b_ptr[idx];
}

#define LT -1
#define GT 1
#define EQ 0

typedef int (*Sorter)(const size_t a, const size_t b,
                      const Instance *restrict const inst);

int sort_max_price(const size_t a, const size_t b,
                   const Instance *restrict const inst) {

    double diff = inst->p_ptr[a] - inst->p_ptr[b];
    if (diff < 0)
        return LT;
    if (diff > 0)
        return GT;
    return EQ;
}

int sort_min_price(const size_t a, const size_t b,
                   const Instance *restrict const inst) {
    return -sort_max_price(a, b, inst);
}

int sort_min_weight(const size_t a, const size_t b,
                    const Instance *restrict const inst) {
    double wa = inst->a_ptr[a] + inst->rho * sqrt(inst->b_ptr[a]);
    double wb = inst->a_ptr[b] + inst->rho * sqrt(inst->b_ptr[b]);
    double diff = wa - wb;
    if (diff < 0)
        return GT;
    if (diff > 0)
        return LT;
    return EQ;
}

int sort_max_weight(const size_t a, const size_t b,
                    const Instance *restrict const inst) {
    return -sort_min_weight(a, b, inst);
}

int sort_max_ratio(const size_t a, const size_t b,
                   const Instance *restrict const inst) {

    double wa = inst->a_ptr[a] + inst->rho * sqrt(inst->b_ptr[a]);
    double ra = inst->p_ptr[a] / wa;
    double wb = inst->a_ptr[b] + inst->rho * sqrt(inst->b_ptr[b]);
    double rb = inst->p_ptr[b] / wb;

    double diff = ra - rb;
    if (diff < 0)
        return LT;
    if (diff > 0)
        return GT;
    return EQ;
}

int sort_min_ratio(const size_t a, const size_t b,
                   const Instance *restrict const inst) {
    return -sort_max_ratio(a, b, inst);
}

int sort_random(const size_t a, const size_t b,
                const Instance *restrict const inst) {
    (void)a;
    (void)b;
    (void)inst;

    int vals[] = {LT, GT, EQ};
    int idx = int_in_range(0, sizeof(vals) / sizeof(int));
    return vals[idx];
}

typedef struct {
    size_t *buff;
    size_t len;
    size_t cap;
} Queue;

int queue_init(Queue *restrict queue, size_t cap) {

    size_t *buff = malloc(sizeof(size_t) * cap);
    if (buff == NULL)
        return 1;
    queue->buff = buff;

    queue->len = 0;
    queue->cap = cap;

    return 0;
}

void queue_swap(Queue *restrict const queue, size_t a, size_t b) {
    int tmp = queue->buff[a];
    queue->buff[a] = queue->buff[b];
    queue->buff[b] = tmp;
}

void heapify_up(Queue *restrict const queue,
                const Instance *restrict const inst, Sorter sorter) {

    bool run = true;
    size_t curr = queue->len;
    while (run) {
        size_t parent = curr / 2;
        int order = sorter(queue->buff[curr], queue->buff[parent], inst);
        if (order == GT) {
            queue_swap(queue, curr, parent);
            curr = parent;
        } else {
            run = false;
        }

        if (curr == 0)
            run = false;
    }
}

void queue_push(Queue *restrict const queue,
                const Instance *restrict const inst, Sorter sorter,
                const size_t item) {
    queue->buff[queue->len] = item;
    heapify_up(queue, inst, sorter);
    queue->len++;
}

bool queue_peak(const Queue *restrict const queue) { return queue->len > 0; }

void heapify_down(Queue *restrict const queue,
                  const Instance *restrict const inst, Sorter sorter) {

    bool run = true;
    size_t curr = 0;
    while (run) {
        size_t left = 2 * curr + 1;
        size_t right = 2 * curr + 2;
        if (left < queue->len &&
            sorter(queue->buff[curr], queue->buff[left], inst) == LT) {
            queue_swap(queue, curr, left);
            curr = left;
        } else if (right < queue->len &&
                   sorter(queue->buff[curr], queue->buff[right], inst) == LT) {
            queue_swap(queue, curr, left);
            curr = right;
        } else {
            run = false;
        }
    }
}

size_t queue_pop(Queue *restrict const queue,
                 const Instance *restrict const inst, Sorter sorter) {
    size_t output = queue->buff[0];
    queue->len--;
    queue->buff[0] = queue->buff[queue->len];
    heapify_down(queue, inst, sorter);
    return output;
}

void queue_drain(Queue *restrict const queue) { queue->len = 0; }

void queue_free(Queue *restrict const queue) { free(queue->buff); }

typedef void (*Remove)(const TabuTenure *restrict, Solution *restrict const,
                       const Instance *restrict const, const int);

void remove_item(const TabuTenure *restrict tt, Solution *restrict const sol,
                 const Instance *restrict const inst, const int idx) {
    sol->items[idx] = false;
    tt->block_move[idx] = tt->block_size;
    sol->a_weight -= inst->a_ptr[idx];
    sol->b_weight -= inst->b_ptr[idx];
    sol->obj_val -= inst->p_ptr[idx];
}

void add_item(const TabuTenure *restrict tt, Solution *restrict const sol,
              const Instance *restrict const inst, int idx) {

    sol->items[idx] = true;
    tt->block_move[idx] = tt->block_size;
    sol->a_weight += inst->a_ptr[idx];
    sol->b_weight += inst->b_ptr[idx];
    sol->obj_val += inst->p_ptr[idx];
}

bool can_insert_item(const TabuTenure *restrict const tt,
                     const Solution *restrict const sol,
                     const Instance *restrict const inst, const size_t idx) {
    double new_weight = (sol->a_weight + inst->a_ptr[idx]) +
                        inst->rho * sqrt(sol->b_weight + inst->b_ptr[idx]);

    return tt_is_free(tt, idx) && !sol->items[idx] &&
           (new_weight <= inst->capacity);
}

typedef struct {
    TabuTenure *tt;
    Solution *sol;
    Queue *queue;
    Sorter sorter;
} OperatorInput;

void remove_items(OperatorInput input, const Instance *restrict const inst,
                  int count) {

    for (size_t i = 0; i < inst->n_items; i++) {
        if (tt_is_free(input.tt, i) && input.sol->items[i]) {
            queue_push(input.queue, inst, input.sorter, i);
        }
    }

    while (count-- && queue_peak(input.queue)) {
        int curr = queue_pop(input.queue, inst, input.sorter);
        remove_item(input.tt, input.sol, inst, curr);
    }
    queue_drain(input.queue);
}

void insert_items(OperatorInput input, const Instance *restrict const inst) {
    for (size_t i = 0; i < inst->n_items; i++) {
        if (can_insert_item(input.tt, input.sol, inst, i)) {
            queue_push(input.queue, inst, input.sorter, i);
        }
    }

    while (queue_peak(input.queue)) {
        int curr = queue_pop(input.queue, inst, input.sorter);
        if (can_insert_item(input.tt, input.sol, inst, curr)) {
            add_item(input.tt, input.sol, inst, curr);
        }
    }

    queue_drain(input.queue);
}

Sorter get_remove_sorter(void) {

    Sorter sorters[] = {sort_min_price, sort_max_weight, sort_min_ratio,
                        sort_random};
    int size = sizeof(sorters) / sizeof(Sorter);
    int idx = int_in_range(0, size);
    return sorters[idx];
}
Sorter get_insert_sorter(void) {
    Sorter sorters[] = {sort_max_price, sort_min_weight, sort_max_ratio};
    int size = sizeof(sorters) / sizeof(Sorter);
    int idx = int_in_range(0, size);
    return sorters[idx];
}

void apply_move(Queue *restrict const queue,
                const Instance *restrict const inst, Solution *sol,
                TabuTenure *restrict tt) {

    int min_update = max(inst->n_items / MAX_UPDATE_SCALE, 1);
    int max_update = max(inst->n_items / MIN_UPDATE_SCALE, 1);
    int count = int_in_range(min_update, max_update);

    OperatorInput op_input = {0};
    op_input.tt = tt;
    op_input.queue = queue;
    op_input.sol = sol;

    op_input.sorter = get_remove_sorter();
    remove_items(op_input, inst, count);

    op_input.sorter = get_insert_sorter();
    insert_items(op_input, inst);
}

void clone_solution(Solution *src, Solution *dst) {
    dst->a_weight = src->a_weight;
    dst->b_weight = src->b_weight;
    dst->obj_val = src->obj_val;
    for (size_t i = 0; i < src->n_items; i++) {
        dst->items[i] = src->items[i];
    }
}

void update_best(Solution *curr, Solution *best) {

    if (curr->obj_val > best->obj_val) {
        clone_solution(curr, best);
    }
}

typedef struct {
    Solution *curr;
    Solution *best;
    TabuTenure *tt;
    Queue *queue;
} HeurEngine;

void tabu_search(HeurEngine *const heur, const Instance *restrict const inst) {

    size_t n_iters = inst->n_items;

    double prev_sol = -1;
    size_t same_sol_count = 0;
    size_t max_same_sol_iters = MAX_SAME_SOL;

    for (size_t i = 0; i < n_iters; i++) {
        for (size_t j = 0; j < n_iters; j++) {
            apply_move(heur->queue, inst, heur->curr, heur->tt);
            update_best(heur->curr, heur->best);
            tt_update(heur->tt);
        }
        double curr = heur->best->obj_val;
        if (prev_sol != curr) {
            same_sol_count = 0;
            prev_sol = curr;

        } else {
            same_sol_count++;
        }

        if (same_sol_count == max_same_sol_iters)
            break;

        tt_reset(heur->tt);
    }
}

void find_initial_solution(Solution *restrict const curr_sol,
                           Solution *restrict const best_sol,
                           const Instance *restrict const inst) {
    solution_init(curr_sol, inst->n_items);
    initial_solution(curr_sol, inst);
    solution_init_from(curr_sol, best_sol);
}

double heuristic_solution(const Instance* const RESTRICT inst, int* ts_sol_items) {

    srand(42);

    Queue queue = { 0 };
    queue_init(&queue, inst->n_items);

    Solution curr_sol = { 0 };
    Solution best_sol = { 0 };
    find_initial_solution(&curr_sol, &best_sol, inst);
    double init_sol = curr_sol.obj_val;

    TabuTenure tabu_tenure = { 0 };
    tabu_tenure_init(&tabu_tenure, inst->n_items);

    HeurEngine engine = { 0 };
    engine.curr = &curr_sol;
    engine.best = &best_sol;
    engine.queue = &queue;
    engine.tt = &tabu_tenure;

    tabu_search(&engine, inst);

    solution_free(&curr_sol);
    for (int i = 0; i < inst->n_items; ++i)
        ts_sol_items[i] = best_sol.items[i];
    solution_free(&best_sol);
    tt_free(&tabu_tenure);
    queue_free(&queue);

    return (int)best_sol.obj_val;
}
