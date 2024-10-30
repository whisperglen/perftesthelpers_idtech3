
#ifndef _CSV_H
#define _CSV_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

static FILE* local_f_csv = NULL;

void csv_open(const char* name);
void csv_close();

void csv_put_string(const char* str);
void csv_put_float(double val);


#ifdef __cplusplus
}
#endif
#endif