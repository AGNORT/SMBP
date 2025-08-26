#ifndef CLI_PARSE_H
#define CLI_PARSE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

typedef struct {
    const char *prog_name;
	const char *input_file;		//the input file
	const char *config_file;	//redundant parameter
	const char *output_file;	//the output file
	const char *problem_type;   //the type of the problem to be solved, including "SMKP" and "SMBP"
	const char *method;         //the method to solve the problem, default is DP or BPC, specifed as Gurobi to use the gurobi solver
	const char *algoCombination; //the combination of algorithms for SMBP, including ULA-VD, ULA-FD, LLA-VD, LLA-FD, default use ULA-VD
	bool SR3;					//To indicate if use SR3 inequalities, default use SR3
	bool DIs;					//To indicate if use DIs inequalities, default not use DIs
    bool direct_rho;            //need this to be true if solve the instances from Xu et al.(2023)
} Args;

int parse_args(Args* args, int argc, char** argv);

#ifdef __cplusplus
}
#endif

#endif
