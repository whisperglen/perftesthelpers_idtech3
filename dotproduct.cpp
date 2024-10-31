#include <math.h>
#include <float.h>
#include <cassert>
#include <intrin.h>
#include <stdio.h>

#include "bdmpx.h"
#include "timing.h"
#include "platform.h"
#include "csv.h"

typedef float vec_t;
typedef vec_t vec2_t[2];
typedef vec_t vec3_t[3];
typedef vec_t vec4_t[4];
typedef vec_t vec5_t[5];

#if idneon
static ID_INLINE float DotProduct_neon(const vec3_t a, const vec3_t b)
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

#if idsse
#define zero_vec 0

static ID_INLINE float DotProduct_sse(const vec4_t a, const vec4_t b)
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

static ID_INLINE float DotProduct_ssev2(const vec4_t a, const vec4_t b)
{
	float ret;
	__m128 aVec, bVec, zero;
	__m128i clr;
#if zero_vec
	zero = _mm_setzero_ps();
#else
#define zero aVec
#endif
  clr = _mm_setr_epi8(0,1,2,3,4,5,6,7,8,9,10,11,(char)128,(char)128,(char)128,(char)128);
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

static ID_INLINE float DotProduct_ssev3(const vec4_t a, const vec4_t b)
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
#endif

static ID_INLINE float DotProduct(const vec4_t x, const vec4_t y)
{
  return ((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2]);
}

#define FUNCTION_INLINED(NAME) \
void finline_##NAME(int reps, int inputsz, const vec4_t* inx, const vec4_t* iny, float *outputs, double *elapsed) \
{ \
	Timer timer; \
	for (int r = 0; r < reps; r++) \
		for (int i = 0; i < inputsz; i++) \
		{ \
			outputs[i] = NAME(inx[i], iny[i]); \
		} \
	*elapsed = timer.elapsed_ms(); \
}
typedef void (*dotp_fn_inlined)(int reps, int inputsz, const vec4_t* inx, const vec4_t* iny, float* outputs, double* elapsed);
typedef float (*dotp_fn)(const vec4_t x, const vec4_t y);

struct function_data
{
	void* fp;
	const char* fname;
};

FUNCTION_INLINED(DotProduct);
#if idsse
FUNCTION_INLINED(DotProduct_sse);
FUNCTION_INLINED(DotProduct_ssev2);
FUNCTION_INLINED(DotProduct_ssev3);
#endif
#if idneon
FUNCTION_INLINED(DotProduct_neon);
#endif

static struct function_data data_fn_inlined[] =
{
	{ finline_DotProduct,        "dotp_scalar" },
#if idsse
	{ finline_DotProduct_sse,    "   dotp_sse" },
	{ finline_DotProduct_ssev2,  "dotp_sse_v2" },
	{ finline_DotProduct_ssev3,  "dotp_sse_v3" },
#endif
#if idneon
	{ finline_DotProduct_neon,   "  dotp_neon" },
#endif
	{ 0, 0 }
};


#define ARRAY_SIZE(X) (sizeof(X)/sizeof(X[0]))

//theres about 7.300.000 samples in the datafile
//there is alot of data here, maybe we care about the cache sizes as well ;-)
//#define TEST_SIZE (128*1000)
#define TEST_SIZE (8*1000*1000)

static vec4_t ina[TEST_SIZE];
static vec4_t inb[TEST_SIZE];
static float output[(ARRAY_SIZE(data_fn_inlined) - 1) * TEST_SIZE];

static int VectorCompare( const vec3_t v1, const vec3_t v2 ) {
	if (v1[0] != v2[0] || v1[1] != v2[1] || v1[2] != v2[2]) {
		return 0;
	}			
	return 1;
}

#define REPETITIONS 1
void maintest_dotproduct(void)
{
  int i,j,k,sz = 0;
  int tested = 0;
  int ercd = 0;

  double elapsed[ARRAY_SIZE(data_fn_inlined) - 1];

  memset(output, 0, sizeof(output));

  ercd = csv_open("./results.csv");
  if (ercd != 0)
  {
	  printf("Could not open the csv file.\n\n");
  }

  //bdmpx_set_option(BDMPX_OPTION_PRECACHE_FILE_FOR_READ, 1);

  ercd = bdmpx_create(NULL, "bdmpx_dotproduct.bin", BDMPX_OP_READ);
  if (ercd != 0)
  {
	  printf("Could not open the test input file. Aborting.\n");
	  return;
  }

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
  printf("Test inputs %d\n", j);

  Timer timer;
  float* out;

  for (tested = 0, out = output; data_fn_inlined[tested].fp != 0; tested++, out += TEST_SIZE)
  {
	  dotp_fn_inlined callme = (dotp_fn_inlined)data_fn_inlined[tested].fp;
	  const char* myinfo = data_fn_inlined[tested].fname;

	  callme(REPETITIONS, j, ina, inb, out, &elapsed[tested]);
	  printf("%d %s: %4.4f\n", tested, myinfo, elapsed[tested]);
	  csv_put_float(elapsed[tested]);
  }

  csv_put_string(",");
  for (k = 1; k < tested; k++)
	  csv_put_float(elapsed[k] / elapsed[0]);

  csv_put_string("\n");
  csv_close();

  for (k = 1; k < tested; k++)
  {
	  const char* myinfo = data_fn_inlined[k].fname;
	  for (i = 0; i < j; i++)
	  {
		  if (output[i] != output[k*TEST_SIZE + i])
			  printf("%d %s: %4.4f %4.4f\n", i, myinfo, output[i], output[k * TEST_SIZE + i]);
	  }
  }
}
