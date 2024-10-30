#include <math.h>
#include <float.h>
#include <cassert>
#include <intrin.h>
#include <stdio.h>

#include "bdmpx.h"
#include "timing.h"
#include "platform.h"

typedef float vec_t;
typedef vec_t vec2_t[2];
typedef vec_t vec3_t[3];
typedef vec_t vec4_t[4];
typedef vec_t vec5_t[5];

#if 0
static float DotProduct_neon(const vec3_t a, const vec3_t b)
{
	float ret = 0;
	float32x2_t aVec, bVec;
	//uint64x1_t iVecS;
	aVec = vld1_f32(a);
	bVec = vld1_f32(b);
	aVec = vmul_f32(aVec, bVec);
	//iVec = vshl_n_u64(aVec, 32);
	//aVec = vadd_f32(aVec, iVec);
	aVec = vpadd_f32(aVec, aVec);
	vst1_lane_f32(&ret, aVec, 0);
	return ret + a[2] * b[2];
}
#endif

#define zero_vec 0

static float DotProduct_sse(const vec4_t a, const vec4_t b)
{
	float ret;
	__m128 aVec, bVec, zero, clr;
#if zero_vec
	zero = _mm_setzero_ps();
#else
#define zero aVec
#endif
	clr = _mm_setr_ps(1., 1., 1., 0.);
	aVec = _mm_loadu_ps(a);
	bVec = _mm_loadu_ps(b);
	aVec = _mm_mul_ps(aVec, clr);
	aVec = _mm_mul_ps(aVec, bVec);
	aVec = _mm_hadd_ps(aVec, zero);
	aVec = _mm_hadd_ps(aVec, zero);
	_mm_store_ss(&ret, aVec);
	return ret;
#undef zero
}

static float DotProduct_sse_alt(const vec4_t a, const vec4_t b)
{
	float ret;
	__m128 aVec, bVec, zero;
  __m128i clr;
#if zero_vec
	zero = _mm_setzero_ps();
#else
#define zero aVec
#endif
  clr = _mm_setr_epi8(0,1,2,3,4,5,6,7,8,9,10,11,128,128,128,128);
	//aVec = _mm_castsi128_ps(_mm_lddqu_si128((__m128i*)a));
	//bVec = _mm_castsi128_ps(_mm_lddqu_si128((__m128i*)b));
	aVec = _mm_loadu_ps(a);
	bVec = _mm_loadu_ps(b);
	aVec = _mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(aVec), clr));
	aVec = _mm_mul_ps(aVec, bVec);
	aVec = _mm_hadd_ps(aVec, zero);
	aVec = _mm_hadd_ps(aVec, zero);
	_mm_store_ss(&ret, aVec);
	return ret;
#undef zero
}

static float DotProduct_sse_altalt(const vec4_t a, const vec4_t b)
{
	float ret;
	__m128 aVec, bVec, zero;
  __m128i clr;
	//zero = _mm_setzero_ps();
	aVec = _mm_loadu_ps(a);
	bVec = _mm_loadu_ps(b);
	aVec = _mm_mul_ps(aVec, bVec);
	aVec = _mm_hadd_ps(aVec, aVec);
	_mm_store_ss(&ret, aVec);
	return ret + a[2] * b[2];
}

static float DotProduct(const vec4_t x, const vec4_t y)
{
  return ((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2]);
}


#define ARRAY_SIZE(X) (sizeof(X)/sizeof(X[0]))


#define TEST_SIZE 1000

static vec3_t ina[TEST_SIZE];
static vec3_t inb[TEST_SIZE];
static float output[4][TEST_SIZE];

static int VectorCompare( const vec3_t v1, const vec3_t v2 ) {
	if (v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2]) {
		return 0;
	}			
	return 1;
}

void maintest_dotproduct(void)
{
  int i,j,k,sz = 0;
  float *out;

  double elapsed;

  //bdmpx_set_option(BDMPX_OPTION_PRECACHE_FILE_FOR_READ, 1);

  bdmpx_create(NULL, "bdmpx_dotproduct.bin", BDMPX_OP_READ);

  for(j = 0; j < TEST_SIZE; j++)
  {
    int sza = sizeof(vec3_t) +1;
    int szb = sizeof(vec3_t) +1;
    int ercd = bdmpx_read(NULL, 2, &sza, &ina[j], &szb, inb[j]);

    if(ercd != 2)
      break;

    assert(sza == sizeof(vec3_t));
    assert(szb == sizeof(vec3_t));
  }
  printf("Test num %d\n", j);

  Timer timer;

  out =  &output[0][0];
  timer.reset();
  for(i = 0; i < j; i++)
  {
    out[i] = DotProduct(ina[i],inb[i]);
  }
  elapsed = timer.elapsed_ms();
  printf("DotProduct: %4.4f\n", elapsed);
  
  out =  &output[1][0];
  timer.reset();
  for(i = 0; i < j; i++)
  {
    out[i] = DotProduct_sse(ina[i],inb[i]);
  }
  elapsed = timer.elapsed_ms();
  printf("DotProduct_sse: %4.4f\n", elapsed);
  
  out =  &output[2][0];
  timer.reset();
  for(i = 0; i < j; i++)
  {
    out[i] = DotProduct_sse_alt(ina[i],inb[i]);
  }
  elapsed = timer.elapsed_ms();
  printf("DotProduct_sse_alt: %4.4f\n", elapsed);

  out =  &output[3][0];
  timer.reset();
  for(i = 0; i < j; i++)
  {
    out[i] = DotProduct_sse_altalt(ina[i],inb[i]);
  }
  elapsed = timer.elapsed_ms();
  printf("DotProduct_sse_alt1: %4.4f\n", elapsed);
#if 0
  timer.reset();
  for(i = 0; i < j; i++)
  {
    output[3][i] = Q_rsqrt_sse_x(input[i]);
  }
  elapsed = timer.elapsed_ms();
  printf("Q_rsqrt_sse_x: %4.4f\n", elapsed);
#endif

  for(i = 0; i < j; i++)
  {
    if(output[0][i] != output[1][i])
      printf("1%d %4.4f %4.4f\n",i,output[0][i],output[1][i]);

    if(output[0][i] != output[2][i])
      printf("2%d %4.4f %4.4f\n",i,output[0][i],output[2][i]);

    if(output[0][i] != output[3][i])
      printf("2%d %4.4f %4.4f\n",i,output[0][i],output[3][i]);
  }
  
}
