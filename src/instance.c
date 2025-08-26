#define _CRT_SECURE_NO_WARNINGS
#include "instance.h"

#include <stdbool.h>
#include <string.h>
//#include <zlib.h>

#define BUFF_SIZE 4096

typedef struct {
    char* a;
    char* b;
    char* c;
} Triple;

void tokenize_line(Triple* RESTRICT triple, char* line) {

    int idx = 0;
    for (char* str = line, *tok; (tok = strtok(str, " ")) != NULL; str = NULL) {
        switch (idx) {
        case 0:
            triple->a = tok;
            break;
        case 1:
            triple->b = tok;
            break;
        case 2:
            triple->c = tok;
            break;
        }

        idx++;
    }
}

#define alloc_vect(dest)                                                       \
    dest = malloc(sizeof(double) * n_items);                                   \
    if (dest == NULL)                                                          \
        return FAILURE;

int parse_header(Instance* RESTRICT inst, const Triple* RESTRICT triple) {

    size_t n_items = atol(triple->a);
    double capacity = atol(triple->b);
    double rho = atof(triple->c);

    inst->n_items = n_items;
    inst->capacity = capacity;
    inst->rho = rho;

    alloc_vect(inst->a_ptr);
    alloc_vect(inst->b_ptr);
    alloc_vect(inst->p_ptr);
    alloc_vect(inst->p_weight);

    return SUCCESS;
}

typedef struct {
    int line;
} ParseBodyState;

int parse_body(ParseBodyState* state, Instance* RESTRICT inst,
    Triple* RESTRICT triple) {

    inst->p_ptr[state->line] = atof(triple->a);
    inst->a_ptr[state->line] = atof(triple->b);
    inst->b_ptr[state->line] = atof(triple->c);
    inst->p_weight[state->line] = atof(triple->a) / (atof(triple->b) + inst->rho*sqrt(atof(triple->c)));

    state->line++;

    return SUCCESS;
}

#define run_or_fail(f)                                                         \
    if (f)                                                                     \
        return FAILURE;

//int instance_parse(Instance * RESTRICT inst,
//                   const char *const RESTRICT file_path) {
//
//    gzFile file = gzopen(file_path, "rb");
//    if (file == NULL)
//        return FAILURE;
//
//    char *line = malloc(BUFF_SIZE * sizeof(char));
//    if (line == NULL)
//        return FAILURE;
//
//    bool first = true;
//    ParseBodyState state = {0};
//    while (gzgets(file, line, BUFF_SIZE) != NULL) {
//        Triple triple = {0};
//        tokenize_line(&triple, line);
//        if (first) {
//            run_or_fail(parse_header(inst, &triple));
//            first = false;
//        } else {
//            run_or_fail(parse_body(&state, inst, &triple));
//        }
//    }
//
//    free(line);
//    gzclose(file);
//    return SUCCESS;
//}

int instance_parse(Instance* RESTRICT inst, const char* const RESTRICT file_path) {
    FILE* file = fopen(file_path, "r");  // 打开普通文本文件
    if (file == NULL)
        return FAILURE;

    char* line = malloc(BUFF_SIZE * sizeof(char));
    if (line == NULL) {
        fclose(file);
        return FAILURE;
    }

    bool first = true;
    ParseBodyState state = { 0 };

    while (fgets(line, BUFF_SIZE, file) != NULL) {
        Triple triple = { 0 };
        tokenize_line(&triple, line);
        if (first) {
            run_or_fail(parse_header(inst, &triple));
            first = false;
        }
        else {
            run_or_fail(parse_body(&state, inst, &triple));
        }
    }

    free(line);
    fclose(file);
    return SUCCESS;
}

void instance_free(Instance* RESTRICT inst) {
    free(inst->a_ptr);
    free(inst->b_ptr);
    free(inst->p_ptr);
    free(inst->p_weight);
}
