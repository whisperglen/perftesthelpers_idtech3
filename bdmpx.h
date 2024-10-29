
#ifndef BDMPX_H_
#define BDMPX_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BDMPDX_CONFIG * bdmpx_handle;

#define BDMPX_OP_READ 0
#define BDMPX_OP_WRITE 1

int bdmpx_create(bdmpx_handle *ret_hndl, const char *filename, int operation);

int bdmpx_close(bdmpx_handle hndl);

int bdmpx_write(bdmpx_handle hndl, int numparams, ...);

int bdmpx_read(bdmpx_handle hndl, int numparams, ...);

int bdmpx_read_alloc(bdmpx_handle hndl, int numparams, ...);

int bdmpx_rewind(bdmpx_handle hndl);

#define BDMPX_OPTION_PRECACHE_FILE_FOR_READ 1

void bdmpx_set_option(int option, int value);

#ifdef __cplusplus
}
#endif
#endif