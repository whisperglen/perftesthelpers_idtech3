
#include "csv.h"

#include <stdio.h>

//FILE* local_f_csv = NULL;

int csv_open(const char* name)
{
	local_f_csv = fopen(name, "a");
	return ( ! local_f_csv);
}

#define RETRY_ONCE(X) \
do { \
	int retry = 1; \
	while((1 != (X)) && (retry-- > 0)); \
} while(0)

void csv_put_string(const char* str)
{
	if (local_f_csv)
	{
		int sz = strlen(str);
		RETRY_ONCE(fwrite(str, sz, 1, local_f_csv));
		fflush(local_f_csv);
	}
}

void csv_put_float(double val)
{
	if (local_f_csv)
	{
		char str[1024];
		int sz = snprintf(str, sizeof(str), "%4.4f,", val);
		RETRY_ONCE(fwrite(str, sz, 1, local_f_csv));
		fflush(local_f_csv);
	}
}

void csv_close()
{
	if (local_f_csv)
	{
		fclose(local_f_csv);
		local_f_csv = NULL;
	}
}
