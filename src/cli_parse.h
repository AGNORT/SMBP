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
    const char *input_file;
    const char *config_file;
    const char *output_file;
	bool direct_rho;          //To adjust the value of rho
} Args;

int parse_args(Args* args, int argc, char** argv);

#ifdef __cplusplus
}
#endif

#endif
