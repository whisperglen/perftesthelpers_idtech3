#include <math.h>
#include <float.h>
#include <cassert>
#include <stdio.h>
#include <float.h>

#include "platform.h"
#include "bdmpx.h"
#include "timing.h"

#if idsse
#include <intrin.h>
#endif
#if idneon
#include <nenon.h>
#endif

float Q_rsqrt_bitlevel( float number )
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
float Q_rsqrt_sse_euler( float number )
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

	x_f = _mm_mul_ss(tmp0_f, number_f);

	y_f = _mm_rsqrt_ss(number_f);
	x_f = _mm_mul_ss(x_f, y_f);
	x_f = _mm_mul_ss(x_f, y_f);
	x_f = _mm_sub_ss(tmp1_f, x_f);
	x_f = _mm_mul_ss(y_f, x_f);

	_mm_store_ss(&ret, x_f);
	return ret;
	}

float Q_rsqrt_sse( float number )
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
float Q_rsqrt_neon(float number)
{
	float32x2_t number_f;
	float ret;

	if (number == 0.f) return 1.f;

	number_f = vld1_dup_f32(&number);

	number_f = vrsqrte_f32(number_f);

	vst1_lane_f32(&ret, number_f, 0);
	return ret;
}
#endif

float Q_rsqrt_math( float number )
{
	if(number == 0.f) return 1.f/FLT_MIN;
	return 1.f/sqrtf(number);
}

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(X[0]))

float inputs [] =
{
	1.0f,
	250.f,
	0.123f,
	2.8f,
};

float outputs[4][ARRAY_SIZE(inputs)];

#define TEST_SIZE 3 * 1000 * 1000

float input[TEST_SIZE];
float output[4][TEST_SIZE];

#define REPETITIONS 1
void maintest_rsqrt(void)
{
	int i,j,k,sz = 0;
	int tests = 0;

	double elapsed;

	memset(output, 0, sizeof(output));

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
		{ //let's skip over zeroes since the precision comparison gives bogus results
			j++;
		}
	}
	printf("Test num %d\n", j);

	Timer timer;

	timer.reset();
	for(i = 0; i < j; i++)
	{
		output[0][i] = Q_rsqrt_math(input[i]);
	}
	elapsed = timer.elapsed_ms();
	tests++;
	printf("0 Q_rsqrt_math: %4.4f\n", elapsed);

	timer.reset();
	for(i = 0; i < j; i++)
	{
		output[1][i] = Q_rsqrt_bitlevel(input[i]);
	}
	elapsed = timer.elapsed_ms();
	tests++;
	printf("1 Q_rsqrt_bitlevel: %4.4f\n", elapsed);

#if idsse
	timer.reset();
	for(i = 0; i < j; i++)
	{
		output[2][i] = Q_rsqrt_sse_euler(input[i]);
	}
	elapsed = timer.elapsed_ms();
	tests++;
	printf("2 Q_rsqrt_sse_euler: %4.4f\n", elapsed);

	timer.reset();
	for(i = 0; i < j; i++)
	{
		output[3][i] = Q_rsqrt_sse(input[i]);
	}
	elapsed = timer.elapsed_ms();
	tests++;
	printf("3 Q_rsqrt_sse: %4.4f\n", elapsed);
#endif

#if idneon
	timer.reset();
	for(i = 0; i < j; i++)
	{
		output[2][i] = Q_rsqrt_neon(input[i]);
	}
	elapsed = timer.elapsed_ms();
	printf("2 Q_rsqrt_neon: %4.4f\n", elapsed);
#endif

	for(k = 1; k < tests; k++)
	{
		double min = 100.0, max = 0.0, avg = 0.0;
		double diff;

		for(i = 0; i < j; i++)
		{
			diff = fabs((double)output[0][i] - output[k][i]);
			//if(diff > 0.01)
			//printf("step %d %d %g %g %g %g\n",k,i,input[i],output[0][i],output[1][i],output[2][i]);
			if(diff < min) min = diff;
			if(diff > max) max = diff;
			avg += diff;
		}
		printf("diff %d %g %.9f %.9f\n", k, min, max, avg / j);
	}
}
