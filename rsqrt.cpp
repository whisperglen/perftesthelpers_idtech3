#include <math.h>
#include <float.h>
#include <cassert>
#include <stdio.h>
#include <float.h>

#include "platform.h"
#include "bdmpx.h"
#include "timing.h"
#include "csv.h"

#if idsse
#include <intrin.h>
#endif
#if idneon
#include <nenon.h>
#endif

ID_INLINE float Q_rsqrt_q3( float number )
{
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	if (number == 0.f) return 1.f / FLT_MIN;

	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;						// evil floating point bit level hacking
	i  = 0x5f3759df - ( i >> 1 );               // what the fuck?
	y  = * ( float * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
	//y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

	//assert( !_isnan((double)y) ); // bk010122 - FPE?

	return y;
}

#if idsse
ID_INLINE float Q_rsqrt_sse_precise( float number )
{
	__m128 x_f, y_f;
	__m128 number_f;
	__m128 tmp0_f, tmp1_f;
	float ret;

	if (number == 0.f) return 1.f / FLT_MIN;

	//number_f = _mm_load_ss(&number);
	number_f = _mm_set_ss(number);
	tmp0_f = _mm_set_ss(0.5f);
	tmp1_f = _mm_set_ss(1.5f);

	y_f = _mm_rsqrt_ss(number_f);

	x_f = _mm_mul_ss(tmp0_f, number_f);
	x_f = _mm_mul_ss(x_f, y_f);
	x_f = _mm_mul_ss(x_f, y_f);
	x_f = _mm_sub_ss(tmp1_f, x_f);
	x_f = _mm_mul_ss(x_f, y_f);

	_mm_store_ss(&ret, x_f);
	return ret;
}

ID_INLINE float Q_rsqrt_sse( float number )
{
	__m128 number_f;
	float ret;

	if (number == 0.f) return 1.f / FLT_MIN;

	//number_f = _mm_load_ss(&number);
	number_f = _mm_set_ss(number);

	number_f = _mm_rsqrt_ss(number_f);

	_mm_store_ss(&ret, number_f);
	return ret;
}
#endif

#if idneon
ID_INLINE float Q_rsqrt_neon(float number)
{
	float32x2_t number_f;
	float ret;

	if (number == 0.f) return 1.f / FLT_MIN;

	number_f = vld1_dup_f32(&number);

	number_f = vrsqrte_f32(number_f);

	vst1_lane_f32(&ret, number_f, 0);
	return ret;
}
#endif

ID_INLINE float Q_rsqrt_math( float number )
{
	if (number == 0.f) return 1.f / FLT_MIN;

	return 1.f/sqrtf(number);
}

#define FUNCTION_INLINED(NAME) \
void finline_##NAME(int reps, int inputsz, float *inputs, float *outputs, double *elapsed) \
{ \
	Timer timer; \
	for (int r = 0; r < reps; r++) \
		for (int i = 0; i < inputsz; i++) \
		{ \
			outputs[i] = NAME(inputs[i]); \
		} \
	*elapsed = timer.elapsed_ms(); \
}
typedef void (*rsqrt_fn_inlined)(int reps, int inputsz, float* inputs, float* outputs, double* elapsed);
typedef float (*rsqrt_fn)(float number);

struct function_data
{
	void* fp;
	const char* fname;
};

static struct function_data data_fn[] =
{
	{ Q_rsqrt_math,        "       rsqrt_math" },
	{ Q_rsqrt_q3,          "         rsqrt_q3" },
#if idsse
	{ Q_rsqrt_sse_precise, "rsqrt_sse_precise" },
	{ Q_rsqrt_sse,         "        rsqrt_sse" },
#endif
#if idneon
	{ Q_rsqrt_neon,        "       rsqrt_neon" },
#endif
	{ 0, 0 }
};

FUNCTION_INLINED(Q_rsqrt_math);
FUNCTION_INLINED(Q_rsqrt_q3);
#if idsse
FUNCTION_INLINED(Q_rsqrt_sse_precise);
FUNCTION_INLINED(Q_rsqrt_sse);
#endif
#if idneon
FUNCTION_INLINED(Q_rsqrt_neon);
#endif

static struct function_data data_fn_inlined[] =
{
	{ finline_Q_rsqrt_math,        "       inl_rsqrt_math" },
	{ finline_Q_rsqrt_q3,          "         inl_rsqrt_q3" },
#if idsse
	{ finline_Q_rsqrt_sse,         "        inl_rsqrt_sse" },
	{ finline_Q_rsqrt_sse_precise, "inl_rsqrt_sse_precise" },
#endif
#if idneon
	{ finline_Q_rsqrt_neon,        "       inl_rsqrt_neon" },
#endif
	{ 0, 0 }
};

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(X[0]))

float shorttest_in [] =
{
	1.0f,
	250.f,
	0.123f,
	2.8f,
};

float shorttest_out[4][ARRAY_SIZE(shorttest_in)];

#define TEST_SIZE 3 * 1000 * 1000

float input[TEST_SIZE];
float output[4 * TEST_SIZE];

#define REPETITIONS 10
void maintest_rsqrt(void)
{
	int i,j,k,sz = 0;
	int tested = 0;

	memset(output, 0, sizeof(output));

	csv_open("./results.csv");

	//bdmpx_set_option(BDMPX_OPTION_PRECACHE_FILE_FOR_READ, 1);

	int ercd = bdmpx_create(NULL, "bdmpx_rsqrt.bin", BDMPX_OP_READ);
	if (ercd != 0)
	{
		printf("Could not open the test input file. Aborting.\n");
		return;
	}

	for(j = 0; j < TEST_SIZE; )
	{
		if(1 != bdmpx_read(NULL, 1, NULL, &input[j]))
		{
			break;
		}
		if (input[j] != 0.f)
		{
			//let's skip over zeroes since the precision comparison gives bogus results
			j++;
		}
	}
	printf("Test inputs %d\n", j);

	Timer timer;
	float* out;

	printf("\n>Force noinline results:\n");
	for (tested = 0, out = output; data_fn[tested].fp != 0; tested++, out += TEST_SIZE)
	{
		rsqrt_fn callme = (rsqrt_fn)data_fn[tested].fp;
		const char* myinfo = data_fn[tested].fname;
		double elapsed;

		timer.reset();
		for (int r = 0; r < REPETITIONS; r++)
			for (i = 0; i < j; i++)
			{
				out[i] = callme(input[i]);
			}
		elapsed = timer.elapsed_ms();
		printf("%d %s: %4.4f\n", tested, myinfo, elapsed);
		//csv_put_float(elapsed);
	}

	printf("\n>Inline results:\n");
	for (tested = 0, out = output; data_fn_inlined[tested].fp != 0; tested++, out += TEST_SIZE)
	{
		rsqrt_fn_inlined callme = (rsqrt_fn_inlined)data_fn_inlined[tested].fp;
		const char* myinfo = data_fn_inlined[tested].fname;
		double elapsed;

		callme(REPETITIONS, j, input, out, &elapsed);
		printf("%d %s: %4.4f\n", tested, myinfo, elapsed);
		csv_put_float(elapsed);
	}

	csv_put_string("\n");
	csv_close();

	printf("\ndiff: num min max avg\n");
	for(k = 1; k < tested; k++)
	{
		double min = 100.0, max = 0.0, avg = 0.0;
		double diff;

		for(i = 0; i < j; i++)
		{
			diff = fabs((double)output[i] - output[k*TEST_SIZE+i]);
			//if(diff > 0.01)
			//printf("step %d %d %g %g %g %g\n",k,i,input[i],output[0][i],output[1][i],output[2][i]);
			if(diff < min) min = diff;
			if(diff > max) max = diff;
			avg += diff;
		}
		printf("diff: %d %g %.9f %.9f\n", k, min, max, avg / j);
	}
}
