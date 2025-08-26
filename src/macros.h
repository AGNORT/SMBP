#ifndef MACROS_H
#define MACROS_H
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#define array_len(a) sizeof(a) / sizeof(a[0])

#define out_file stdout

#define run_or_fail(f)                                                         \
    if (f)                                                                     \
        return FAILURE;

#define todo(msg)                                                              \
    fprintf(out_file, "missing implementation '%s' in %s line %d", msg,        \
            __FILE__, __LINE__);                                               \
    abort();

#define SUCCESS 0
#define FAILURE 1
#define RUN_OR_FAIL(f)                                                         \
    if (f)                                                                     \
        return FAILURE;
#define MEM_OR_FAIL(f)                                                         \
    if (f == NULL)                                                             \
        return NULL;
#define FORMAT_SPECIFIER(type)                                                 \
    _Generic((type),                                                           \
        char: "%c",                                                            \
        signed char: "%hhd",                                                   \
        unsigned char: "%hhu",                                                 \
        signed short: "%hd",                                                   \
        unsigned short: "%hu",                                                 \
        signed int: "%d",                                                      \
        unsigned int: "%u",                                                    \
        long int: "%ld",                                                       \
        unsigned long int: "%lu",                                              \
        long long int: "%lld",                                                 \
        unsigned long long int: "%llu",                                        \
        float: "%f",                                                           \
        double: "%lf",                                                         \
        long double: "%Lf",                                                    \
        char *: "%s",                                                          \
        void *: "%p",                                                          \
        bool: "%d",                                                            \
        default: "%p")

#define debug(a)                                                               \
    fprintf(out_file, "%s ", __FILE__);                                        \
    fprintf(out_file, "%d ", __LINE__);                                        \
    fprintf(out_file, "%s = ", #a);                                            \
    fprintf(out_file, FORMAT_SPECIFIER(a), a);                                 \
    putchar('\n');

#define print(a)                                                               \
    printf(FORMAT_SPECIFIER(a), a);                                            \
    putchar('\n');

#ifdef __cplusplus
}
#endif
#endif
