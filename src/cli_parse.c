#include "cli_parse.h"


void usage(const char *const prog_name) {
    fprintf(stderr, "%s: usage\n", prog_name);
    fprintf(stderr,
            "%s -i <input-file> [-c config_file] [-o output_file] [-d] | -h\n",
            prog_name);
    fprintf(stderr, "\t-i: specify input file [REQUIRED]\n");
    fprintf(stderr, "\t-o: specify output file\n");
    fprintf(stderr, "\t-c: specify config file\n");
    fprintf(stderr, "\t-d: consider probability from instance file as "
                    "robustness coefficent\n");
    fprintf(stderr, "\t-h: show this help and exit\n");
    exit(1);
}

int parse_args(Args *args, int argc, char **argv) {

    if (argc < 2)
        return FAILURE;

    args->prog_name = argv[0];
    args->direct_rho = false;
    args->method = "MyMethod";
	args->SR3 = true;
	args->DIs = false;

    char ch;
    while ((ch = getopt(argc, argv, "i:o:p:m:a:SDr")) != -1) {
        switch (ch) {
        case 'i':
            args->input_file = optarg;
            break;
        case 'o':
            args->output_file = optarg;
            break;
        case 'p':
            args->problem_type = optarg;
			break;
        case 'm':
            args->method = optarg;
			break;
        case 'a':
            args->algoCombination = optarg;
			break;
        case 'S':
            args->SR3 = false;
            break;
        case 'D':
            args->DIs = true;
			break;
        case 'r':
            args->direct_rho = true;
            break;
        case '?':
            usage(argv[0]);
        }
    }

    return SUCCESS;
}
