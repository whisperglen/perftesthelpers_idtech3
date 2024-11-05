
#ifndef _CSV_H
#define _CSV_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

static FILE* local_f_csv = NULL;

int csv_open(const char* name);
void csv_close();

void csv_put_string(const char* str);
void csv_put_float(double val);
void csv_put_int(int val);


#ifdef __cplusplus
}
#endif
#endif