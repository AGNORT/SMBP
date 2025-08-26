#ifndef __INSTANCE_H
#define __INSTANCE_H
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "macros.h"

typedef struct {
	size_t n_items;			// number of items
	double* a_ptr;			// mean weights
	double* b_ptr;          // variance weights    
	double* p_ptr;          // profit
	double* p_weight;       // profit per weight
	double capacity;		// knapsack capacity
	double rho;				// risk parameter
} Instance;

int instance_parse(Instance *RESTRICT, const char* const RESTRICT);
void instance_free(Instance *RESTRICT);

#ifdef __cplusplus
}
#endif
#endif
