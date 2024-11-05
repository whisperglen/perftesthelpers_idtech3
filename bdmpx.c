
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "bdmpx.h"

struct BDMPDX_CONFIG
{
  FILE *f;
  unsigned int f_sz;
  unsigned char *cache;
  unsigned int cache_pos;
};

static struct BDMPDX_CONFIG default_hndl = { NULL };
static int g_op_precache_file_for_read = 0;

#define MAX_NUM_PARAMS 100

#define FREE(X) do { if(X) { free(X); (X)=NULL; } } while(0)

#define RETRY_ONCE(X) \
do { \
  int retry = 1; \
  while((1 != (X)) && (retry-- > 0)); \
} while(0)

int bdmpx_write(bdmpx_handle hndl, int numparams, ...)
{
  int ret = 0;

  FILE *f;

  va_list args;
  int elm_sz[MAX_NUM_PARAMS];
  unsigned char *elm_data[MAX_NUM_PARAMS];

  int i = 0;
  int j;
  int tmp = 0;

  if(hndl == NULL) hndl = &default_hndl;
  if(hndl->f == NULL) bdmpx_create(NULL,NULL,BDMPX_OP_WRITE);
  if(hndl->f == NULL) return -1;
  f = hndl->f;

  va_start(args, numparams);

  if(numparams > MAX_NUM_PARAMS) numparams = MAX_NUM_PARAMS;
  while(numparams > 0)
  {
    elm_sz[i] = va_arg(args, int);
    elm_data[i] = va_arg(args, unsigned char *);

    i++;
    numparams --;
  }

  va_end(args);

  for (j = 0; j < i; j++)
  {
    RETRY_ONCE(fwrite(&elm_sz[j], sizeof(elm_sz[0]), 1, f));
  }
  RETRY_ONCE(fwrite(&tmp, sizeof(elm_sz[0]), 1, f));

  for (j = 0; j < i; j++)
  {
    RETRY_ONCE(fwrite(elm_data[j], elm_sz[j], 1, f));
  }

  fflush(f);

  return ret;
}

static size_t fread_internal(void * dst, size_t esz, size_t cnt, bdmpx_handle local)
{
  size_t ret = 0;
  if(local->cache)
  {
    int to_read = esz * cnt;
    if(to_read + local->cache_pos > local->f_sz)
      to_read = local->f_sz - local->cache_pos;
    if(to_read != 0)
    {
      memcpy(dst,&local->cache[local->cache_pos],to_read);
      local->cache_pos += to_read;
      ret = 1;
    }
  }
  else
    ret = fread(dst,esz,cnt,local->f);

  return ret;
}

static int fseek_internal(bdmpx_handle local, int off, int origin)
{
  int ret = -1;
  if(local->cache)
  {
    int final = local->cache_pos;
    switch(origin)
    {
    case SEEK_CUR:
      final = (int)local->cache_pos + off;
      break;
    case SEEK_END:
      final = (int)local->f_sz + off;
      break;
    case SEEK_SET:
      final = off;
      break;
    }
    if(final < 0) final = 0;
    if((unsigned int)final > local->f_sz) final = local->f_sz;
    ret = local->cache_pos = final;
  }
  else
  {
	  int ercd = fseek(local->f, off, origin);
	  if (ercd == 0) //success
	  {
		  ret = ftell(local->f);
	  }
  }

  return ret;
}

int bdmpx_read(bdmpx_handle hndl, int numparams, ...)
{
  int ret = 0;

  va_list args;
  int elm_sz[MAX_NUM_PARAMS];

  int j;
  int i;
  int mynumparams = numparams > MAX_NUM_PARAMS ? MAX_NUM_PARAMS : numparams;

  if(hndl == NULL) hndl = &default_hndl;
  if(hndl->f == NULL) bdmpx_create(NULL,NULL,BDMPX_OP_READ);
  if(hndl->f == NULL) return -1;

  va_start(args, numparams);

  for (j = 0; /*j < mynumparams*/; j++)
  {
    int tmp;
    size_t res = fread_internal(&tmp, sizeof(tmp), 1, hndl);
    if(j < MAX_NUM_PARAMS)
    {
      elm_sz[j] = tmp;
    }
    if(tmp == 0 || res == 0)
    {
      break;
    }
  }

  for (i = 0; i < j; i++)
  {
    int *usr_sz = 0;
    unsigned char *usr_data = 0;
    int tmp_insz;

    if(i < numparams)
    {
      usr_sz = va_arg(args, int *);
      usr_data = va_arg(args, unsigned char *);
    }

    ret++;
    tmp_insz = elm_sz[i];
    if(usr_sz) { if(tmp_insz < *usr_sz) *usr_sz = tmp_insz; else tmp_insz = *usr_sz; }
    if(usr_data) { fread_internal(usr_data, tmp_insz, 1, hndl); if(tmp_insz < elm_sz[i]) fseek_internal(hndl, elm_sz[i] - tmp_insz, SEEK_CUR); }
    else fseek_internal(hndl, elm_sz[i], SEEK_CUR);
  }

  for( ; i < numparams; i++)
  {
    int *usr_sz = va_arg(args, int *);
    unsigned char *usr_data = va_arg(args, unsigned char *);

    int tmp_insz = 0;
    if(usr_sz) { tmp_insz = *usr_sz; *usr_sz = 0; }
    //if(usr_data && tmp_insz) memset(usr_data, 0, tmp_insz);
  }

  va_end(args);

  return ret;
}

int bdmpx_peek(bdmpx_handle hndl)
{
	int ret = -1;
	int ercd;
	int elm_sz[MAX_NUM_PARAMS];

	int j;
	int i;
	int mynumparams = MAX_NUM_PARAMS;

	if (hndl == NULL) hndl = &default_hndl;
	if (hndl->f == NULL) bdmpx_create(NULL, NULL, BDMPX_OP_READ);
	if (hndl->f == NULL) return -1;

	int rewpos = fseek_internal(hndl, 0, SEEK_CUR);
	if (rewpos >= 0)
	{

		for (j = 0; /*j < mynumparams*/; j++)
		{
			int tmp;
			size_t res = fread_internal(&tmp, sizeof(tmp), 1, hndl);
			if (j < MAX_NUM_PARAMS)
			{
				elm_sz[j] = tmp;
			}
			if (tmp == 0 || res == 0)
			{
				break;
			}
		}

		ercd = fseek_internal(hndl, rewpos, SEEK_SET);
		if (ercd >= 0)
		{
			ret = j;
		}
	}

	return ret;
}

int bdmpx_read_alloc(bdmpx_handle hndl, int numparams, ...)
{
  int ret = 0;

  va_list args;
  int elm_sz[MAX_NUM_PARAMS];

  int j;
  int i;
  int mynumparams = numparams > MAX_NUM_PARAMS ? MAX_NUM_PARAMS : numparams;

  if(hndl == NULL) hndl = &default_hndl;
  if(hndl->f == NULL) bdmpx_create(NULL,NULL,BDMPX_OP_READ);
  if(hndl->f == NULL) return -1;

  va_start(args, numparams);

  for (j = 0; /*j < mynumparams*/; j++)
  {
    int tmp;
    size_t res = fread_internal(&tmp, sizeof(tmp), 1, hndl);
    if(j < MAX_NUM_PARAMS)
    {
      elm_sz[j] = tmp;
    }
    if(tmp == 0 || res == 0)
    {
      break;
    }
  }

  for (i = 0; i < j; i++)
  {
    int *usr_sz = 0;
    void *usr_data = 0;
    int tmp_insz;
    int alloc = 0;

    if(i < numparams)
    {
      usr_sz = va_arg(args, int *);
      usr_data = va_arg(args, unsigned char *);
    }

    ret++;
    tmp_insz = elm_sz[i];
    if(usr_sz)
    {
      if(*usr_sz == 0)
        alloc = 1;
      else
      {
        if(tmp_insz < *usr_sz) *usr_sz = tmp_insz;
        else tmp_insz = *usr_sz;
      }
    }
    else alloc = 1;

    if(usr_data)
    {
      if(alloc)
      {
        unsigned char **dest = (unsigned char**)usr_data;
        unsigned char* alc = (unsigned char*)malloc(tmp_insz);
        if(alc == NULL) goto skip_read;
        *dest = alc;
        usr_data = alc;
      }
      fread_internal(usr_data, tmp_insz, 1, hndl);
      if(tmp_insz < elm_sz[i]) fseek_internal(hndl, elm_sz[i] - tmp_insz, SEEK_CUR);
    }
    else
    {
      skip_read:
        fseek_internal(hndl, elm_sz[i], SEEK_CUR);
    }
  }

  for( ; i < numparams; i++)
  {
    int *usr_sz = va_arg(args, int *);
    unsigned char *usr_data = va_arg(args, unsigned char *);

    int tmp_insz = 0;
    if(usr_sz) { tmp_insz = *usr_sz; *usr_sz = 0; }
    if(usr_data && tmp_insz) memset(usr_data, 0, tmp_insz);
  }

  va_end(args);

  return ret;
}

int bdmpx_rewind(bdmpx_handle hndl)
{
  int ret = 0;

  if(hndl == NULL) hndl = &default_hndl;
  if(hndl->f == NULL) return -1;

  if(hndl->cache)
  {
    hndl->cache_pos = 0;
  }
  else
  {
    fseek(hndl->f,0,SEEK_SET);
  }

  return ret;
}

#define DEFAULT_FILENAME "bdmpx.bin"

int bdmpx_create(bdmpx_handle *ret_hndl, const char *filename, int operation)
{
  int ret = 0;
  FILE *f;

  bdmpx_handle local = NULL;

  if(ret_hndl == NULL)
  {
	if(default_hndl.f != NULL)
	{
		bdmpx_close(&default_hndl);
	}
    local = &default_hndl;
  }
  else
  {
    local = (bdmpx_handle)calloc(1, sizeof(struct BDMPDX_CONFIG));
  }

  if(local == NULL) return -1;

  if(filename == NULL) filename = DEFAULT_FILENAME;

  f = fopen(filename, operation == BDMPX_OP_READ?"rb":"wb");

  if(f)
  {
    local->f = f;
    if(operation == BDMPX_OP_READ)
    {
      fseek(f, 0, SEEK_END);
      local->f_sz = ftell(f);
      fseek(f, 0, SEEK_SET);
      if(g_op_precache_file_for_read)
      {
        local->cache = (unsigned char *)malloc(local->f_sz);
        if(local->cache)
        {
          RETRY_ONCE(fread(local->cache, local->f_sz, 1, local->f));
          fclose(local->f);
        }
      }
    }
  }
  else
  {
    ret = -2;
    if(local != &default_hndl) FREE(local);
  }

  if(ret_hndl != NULL)
    *ret_hndl = local;

  return ret;
}

int bdmpx_close(bdmpx_handle hndl)
{
  int ret = 0;

  bdmpx_handle local = hndl;

  if(local == NULL)
  {
    if(default_hndl.f != NULL)
    {
      local = &default_hndl;
    }
  }

  if(local == NULL) return -1;

  if(local->f) fclose(local->f);
  FREE(local->cache);
  memset(local,0,sizeof(struct BDMPDX_CONFIG));
  if(local != &default_hndl)
    FREE(local);

  return ret;
}

void bdmpx_set_option(int option, int value)
{
  switch(option)
  {
  case BDMPX_OPTION_PRECACHE_FILE_FOR_READ:
    g_op_precache_file_for_read = 1;
    break;
  default:
    break;
  }
}
