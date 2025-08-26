#ifndef COMPACTMOD_H
#define COMPACTMOD_H

#ifdef __cplusplus
extern "C" {
#endif
#include "cli_parse.h"
#include "instance.h"
#ifdef __cplusplus
}
#endif

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include "gurobi_c++.h"


//call the gurobi solver to solve the compact SOCP model of submodular knapsack problem
void SolveCompactKnapsackModel(Args& args);

//call the gurobi solver to solve the linear relaxation of the knapsack problem
void SolveLinearKnapsack(Instance& instance, int itemIdx, double& preWeight, double& dualProfit);


#endif