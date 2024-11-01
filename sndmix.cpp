
#include <math.h>
#include <float.h>
#include <cassert>
#include <stdio.h>

#include "bdmpx.h"
#include "timing.h"
#include "platform.h"
#include "csv.h"
#include "intrin.h"

static int     *snd_p;
static int      snd_linear_count;
static short   *snd_out;

void S_WriteLinearBlastStereo16( void )
{
	int		i;
	int		val;
	int		*src = snd_p;
	short	*dst = snd_out;

	for ( i = 0; i < snd_linear_count; i++, src++, dst++ )
	{
		val = *src>>8;
		if ( val > 32767 )
			*dst = 32767;
		else if ( val < -32768 )
			*dst = -32768;
		else
			*dst = val;
	}
}

#if id386
void S_WriteLinearBlastStereo16_MMX( void )
{
__asm {
	push ebx
	push esi
	push edi
	mov esi,snd_p
	mov edi,snd_out
	mov ebx,snd_linear_count
	test ebx,ebx
	jz	LExit
	mov ecx,esi
	and ecx,63
	jz LMain
	and ecx,3
	jnz LTail
	shr ecx,2
	not ecx
	add ecx,17
LClamp1:
	mov eax,[esi]
	sar eax,8
	cmp eax,32767
	jg	LClampHigh1
	cmp eax,-32768
	jnl LClampDone1
	mov eax,-32768
	jmp LClampDone1
LClampHigh1:
	mov eax,32767
LClampDone1:
	mov [edi],ax
	add esi,4
	add edi,2
	dec ebx
	jz	LExit
	dec ecx
	jnz	LClamp1
LMain:
	mov ecx,ebx
	shr ecx,4
	jz  LTail
	and ebx,15
LAgain:
	movq mm0, qword ptr [esi+ 0]
	movq mm1, qword ptr [esi+ 8]
	movq mm2, qword ptr [esi+16]
	movq mm3, qword ptr [esi+24]
	movq mm4, qword ptr [esi+32]
	movq mm5, qword ptr [esi+40]
	movq mm6, qword ptr [esi+48]
	movq mm7, qword ptr [esi+56]
	psrad mm0,8
	psrad mm1,8
	psrad mm2,8
	psrad mm3,8
	psrad mm4,8
	psrad mm5,8
	psrad mm6,8
	psrad mm7,8
	packssdw mm0, mm1
	packssdw mm2, mm3
	packssdw mm4, mm5
	packssdw mm6, mm7
	movq qword ptr [edi+ 0], mm0
	movq qword ptr [edi+ 8], mm2
	movq qword ptr [edi+16], mm4
	movq qword ptr [edi+24], mm6
	add esi, 64
	add edi, 32
	dec ecx
	jnz LAgain
LTail:
	test ebx, ebx
	jz	LEnd
LClamp2:
	mov eax,[esi]
	sar eax,8
	cmp eax,32767
	jg	LClampHigh2
	cmp eax,-32768
	jnl LClampDone2
	mov eax,-32768
	jmp LClampDone2
LClampHigh2:
	mov eax,32767
LClampDone2:
	mov [edi],ax
	add esi,4
	add edi,2
	dec ebx
	jnz	LClamp2
LEnd:
    emms
LExit:
	pop edi
	pop esi
	pop ebx
} // __asm
}
#endif

#if id386 && idsse
void S_WriteLinearBlastStereo16_SSEasm( void )
{
__asm {
	push ebx
	push esi
	push edi
	mov esi,snd_p
	mov edi,snd_out
	mov ebx,snd_linear_count
	test ebx,ebx
	jz	LExit
	mov ecx,esi
	and ecx,63
	jz LMain
	and ecx,3
	jnz LTail
	shr ecx,2
	not ecx
	add ecx,17
LClamp1:
	mov eax,[esi]
	sar eax,8
	cmp eax,32767
	jg	LClampHigh1
	cmp eax,-32768
	jnl LClampDone1
	mov eax,-32768
	jmp LClampDone1
LClampHigh1:
	mov eax,32767
LClampDone1:
	mov [edi],ax
	add esi,4
	add edi,2
	dec ebx
	jz	LExit
	dec ecx
	jnz	LClamp1
LMain:
	mov ecx,ebx
	shr ecx,5
	jz  LTail
	and ebx,31
LAgain:
	movdqa xmm0, xmmword ptr [esi+0]
	movdqa xmm1, xmmword ptr [esi+16]
	movdqa xmm2, xmmword ptr [esi+32]
	movdqa xmm3, xmmword ptr [esi+48]
	movdqa xmm4, xmmword ptr [esi+64]
	movdqa xmm5, xmmword ptr [esi+80]
	movdqa xmm6, xmmword ptr [esi+96]
	movdqa xmm7, xmmword ptr [esi+112]
	psrad xmm0,8
	psrad xmm1,8
	psrad xmm2,8
	psrad xmm3,8
	psrad xmm4,8
	psrad xmm5,8
	psrad xmm6,8
	psrad xmm7,8
	packssdw xmm0, xmm1
	packssdw xmm2, xmm3
	packssdw xmm4, xmm5
	packssdw xmm6, xmm7
	movdqu xmmword ptr [edi+ 0], xmm0
	movdqu xmmword ptr [edi+16], xmm2
	movdqu xmmword ptr [edi+32], xmm4
	movdqu xmmword ptr [edi+48], xmm6
	add esi, 128
	add edi, 64
	dec ecx
	jnz LAgain
LTail:
	test ebx, ebx
	jz	LEnd
LClamp2:
	mov eax,[esi]
	sar eax,8
	cmp eax,32767
	jg	LClampHigh2
	cmp eax,-32768
	jnl LClampDone2
	mov eax,-32768
	jmp LClampDone2
LClampHigh2:
	mov eax,32767
LClampDone2:
	mov [edi],ax
	add esi,4
	add edi,2
	dec ebx
	jnz	LClamp2
LEnd:
LExit:
	pop edi
	pop esi
	pop ebx
} // __asm
}
#endif

#if idsse

#define SSEALIGNED(X) (((intptr_t)(X) & 15u) == 0)

void S_WriteLinearBlastStereo16_SSE(void)
{
	int		i;
	int		val;
	int		*src = snd_p;
	short	*dst = snd_out;

	__m128i xm0, xm1, xm2, xm3/*, xm4, xm5, xm6, xm7*/;

	int vec_count;

	for (i = 0; !SSEALIGNED(src) && (i < snd_linear_count); i++, src++, dst++)
	{
		val = *src >> 8;
		if (val > 32767)
			*dst = 32767;
		else if (val < -32768)
			*dst = -32768;
		else
			*dst = val;
	}

	vec_count = snd_linear_count - 16;

	for ( ; i <= vec_count; i += 16, src += 16, dst += 16)
	{
		xm0 = _mm_load_si128((__m128i*)src);
		xm1 = _mm_load_si128((__m128i*)src+1);
		xm2 = _mm_load_si128((__m128i*)src+2);
		xm3 = _mm_load_si128((__m128i*)src+3);
		//xm4 = _mm_load_si128((__m128i*)src+4);
		//xm5 = _mm_load_si128((__m128i*)src+5);
		//xm6 = _mm_load_si128((__m128i*)src+6);
		//xm7 = _mm_load_si128((__m128i*)src+7);

		xm0 = _mm_srai_epi32(xm0, 8);
		xm1 = _mm_srai_epi32(xm1, 8);
		xm2 = _mm_srai_epi32(xm2, 8);
		xm3 = _mm_srai_epi32(xm3, 8);
		//xm4 = _mm_srai_epi32(xm4, 8);
		//xm5 = _mm_srai_epi32(xm5, 8);
		//xm6 = _mm_srai_epi32(xm6, 8);
		//xm7 = _mm_srai_epi32(xm7, 8);

		xm0 = _mm_packs_epi32(xm0, xm1);
		xm2 = _mm_packs_epi32(xm2, xm3);/*
		xm4 = _mm_packs_epi32(xm4, xm5);
		xm6 = _mm_packs_epi32(xm6, xm7);*/

		_mm_storeu_si128((__m128i*)dst, xm0);
		_mm_storeu_si128((__m128i*)dst+1, xm2);/*
		_mm_storeu_si128((__m128i*)dst+2, xm4);
		_mm_storeu_si128((__m128i*)dst+3, xm6);*/
	}

	for ( ; i < snd_linear_count; i++, src++, dst++)
	{
		val = *src >> 8;
		if (val > 32767)
			*dst = 32767;
		else if (val < -32768)
			*dst = -32768;
		else
			*dst = val;
	}
}
#undef SSEALIGNED
#endif

#if idneon

#define NEONALIGNED(X) (((intptr_t)(X) & 15u) == 0)

void S_WriteLinearBlastStereo16_neon(void)
{
	int		i;
	int		val;
	int		*src = snd_p;
	short	*dst = snd_out;

	int32x4_t q0, q1, q2, q3;
	int16x4_t d0, d1, d2, d3;
	int16x8_t q4, q5;

	int vec_count;

	for (i = 0; !NEONALIGNED(src) && (i < snd_linear_count); i++, src++, dst++)
	{
		val = *src >> 8;
		if (val > 32767)
			*dst = 32767;
		else if (val < -32768)
			*dst = -32768;
		else
			*dst = val;
	}

	vec_count = snd_linear_count - 16;

	for (; i <= vec_count; i += 16, src += 16, dst += 16)
	{
		q0 = vld1q_s32((const int32_t*)src);
		q1 = vld1q_s32((const int32_t*)src + 4);
		q2 = vld1q_s32((const int32_t*)src + 8);
		q3 = vld1q_s32((const int32_t*)src + 12);

		q0 = vshrq_n_s32(q0, 8);
		q1 = vshrq_n_s32(q1, 8);
		q2 = vshrq_n_s32(q2, 8);
		q3 = vshrq_n_s32(q3, 8);

		d0 = vqmovn_s32(q0);
		d1 = vqmovn_s32(q1);
		d2 = vqmovn_s32(q2);
		d3 = vqmovn_s32(q3);

		q4 = vcombine_s16(d0, d1);
		q5 = vcombine_s16(d2, d3);

		vst1q_s16((int16_t*)dst, q4);
		vst1q_s16((int16_t*)dst + 8, q5);
	}

	for (; i < snd_linear_count; i++, src++, dst++)
	{
		val = *src >> 8;
		if (val > 32767)
			*dst = 32767;
		else if (val < -32768)
			*dst = -32768;
		else
			*dst = val;
	}
}
#undef NEONALIGNED

#endif

#if idneon
void myfunc(const short *samples, int *out)
{
#if 1
	int16x4_t samples_short = vld1_s16(samples); samples +=4;

	  int32x4_t samples_int = vmovl_s16(samples_short);

	  int64x2_t samples_0 = vmovl_s32(vget_low_s32(samples_int));
	int64x2_t samples_0t = vshlq_n_s64(samples_0, 32);
	samples_0 = vaddq_s64(samples_0, samples_0t);

	int64x2_t samples_1 = vmovl_s32(vget_high_s32(samples_int));
	int64x2_t samples_1t = vshlq_n_s64(samples_1, 32);
	samples_1 = vaddq_s64(samples_1, samples_1t);

	  vst1q_s64((int64_t*)out, samples_0);
	  vst1q_s64((int64_t*)out + 2, samples_1);
#elif 1
  int16x4_t samples_short = vld1_s16(samples); samples +=4;
  int32x4_t samples_int = vmovl_s16(samples_short);

  int32x2_t samples_0l = vdup_lane_s32(vget_low_s32(samples_int), 0);
  int32x2_t samples_0r = vdup_lane_s32(vget_low_s32(samples_int), 1);
  int32x2_t samples_1l = vdup_lane_s32(vget_high_s32(samples_int), 0);
  int32x2_t samples_1r = vdup_lane_s32(vget_high_s32(samples_int), 1);

  int32x4_t out_0 = vcombine_s32(samples_0l, samples_0r);
  int32x4_t out_1 = vcombine_s32(samples_1l, samples_1r);

  vst1q_s32(out, out_0);
  vst1q_s32(out + 4, out_1);
#else
  int16x4_t samples_short = vld1_s16(samples); samples +=4;

  int8x8_t swizzle_0 = { 0,1,8,8, 0,1,8,8};
  int8x8_t samples_0 = vtbl1_s8((int8x8_t)samples_short, swizzle_0);
  int8x8_t swizzle_1 = { 2,3,8,8, 2,3,8,8};
  int8x8_t samples_1 = vtbl1_s8((int8x8_t)samples_short, swizzle_1);

  int8x8_t swizzle_2 = { 4,5,8,8, 4,5,8,8};
  int8x8_t samples_2 = vtbl1_s8((int8x8_t)samples_short, swizzle_2);
  int8x8_t swizzle_3 = { 6,7,8,8, 6,7,8,8};
  int8x8_t samples_3 = vtbl1_s8((int8x8_t)samples_short, swizzle_3);

  int32x4_t out_0 = vcombine_s32((int32x2_t)samples_0, (int32x2_t)samples_1);
  int32x4_t out_1 = vcombine_s32((int32x2_t)samples_2, (int32x2_t)samples_3);
  vst1q_s32(out, out_0);
  vst1q_s32(out + 4, out_1);
#endif
}
#endif

typedef float (*sndmix_fn)(void);

struct function_data
{
	void* fp;
	const char* fname;
};

static struct function_data data_fn[] =
{
	{ S_WriteLinearBlastStereo16,        "sndmix_scalar" },
#if id386
	{ S_WriteLinearBlastStereo16_MMX,    "   sndmix_mmx" },
#endif
#if  id386 && idsse
	{ S_WriteLinearBlastStereo16_SSEasm, "sndmix_sseasm" },
#endif
#if idsse
	{ S_WriteLinearBlastStereo16_SSE,    "   sndmix_sse" },
#endif
#if idneon
	{ S_WriteLinearBlastStereo16_neon,   "  sndmix_neon" },
#endif
	{ 0, 0 }
};

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(X[0]))

#define MAX_SAMPLES 100000
static QALIGNA(64)  int     snd_p_d[MAX_SAMPLES]  QALIGNB(64);
static short   snd_out_d[(ARRAY_SIZE(data_fn) -1) * MAX_SAMPLES];

#define TEST_SEGMENTS 900000
#define MAX_ERRORS 10
#define REPETITIONS 1
void maintest_sndmix(void)
{
	//short buff[4]; int buffout[8];
	//buff[0] = 1;
	//buff[1] = 2;
	//buff[2] = 3;
	//buff[3] = 4;
	//myfunc(buff,buffout);
	//printf("%x %x %x %x %x %x %x %x\n",buffout[0],buffout[1],buffout[2],buffout[3],buffout[4],buffout[5],buffout[6],buffout[7]);


#if 1
  int i, j, k;
  int tested = 0;
  int err = 0;
  Timer timer[ARRAY_SIZE(data_fn) - 1];
  Cputime cpu[ARRAY_SIZE(data_fn) - 1];
  double elapsed[ARRAY_SIZE(data_fn) - 1];
  double cputime[ARRAY_SIZE(data_fn) - 1];

  memset(snd_p_d, 0, sizeof(snd_p_d));

  err = csv_open("./results_sndmix.csv");
  if (err != 0)
  {
	  printf("Could not open the csv file.\n\n");
  }

  err = bdmpx_create(NULL, "bdmpx_sndmix.bin", BDMPX_OP_READ);
  if (err != 0)
  {
	  printf("Could not open the test input file. Aborting.\n");
	  return;
  }

  for(j = 0; j < TEST_SEGMENTS && err < MAX_ERRORS ; j++)
  {
    int szlc = sizeof(int)+1;
    int szsrc = sizeof(snd_p_d)+1;
	int samples = 0;
    //bdmpx_write(NULL, 2, sizeof(int), &snd_linear_count, sizeof(int)*snd_linear_count, src);
    int ercd = bdmpx_read(NULL, 2, &szlc, &samples, &szsrc, snd_p_d);

    if(ercd != 2)
      break;
    assert(szlc == sizeof(int));
    assert(szsrc == sizeof(int) * samples);

	{
		//dry run
		snd_p = snd_p_d;
		snd_out = &snd_out_d[0];
		snd_linear_count = samples;
		S_WriteLinearBlastStereo16();
	}

	for (tested = 0; data_fn[tested].fp != 0; tested++)
	{
		sndmix_fn callme = (sndmix_fn)data_fn[tested].fp;
		const char* myinfo = data_fn[tested].fname;

		snd_p = snd_p_d;
		snd_out = &snd_out_d[tested * MAX_SAMPLES];
		snd_linear_count = samples;

		cpu[tested].reset();
		timer[tested].reset();
		for (int r = 0; r < REPETITIONS; r++)
			callme();

		elapsed[tested] = timer[tested].accum();
		cputime[tested] = cpu[tested].accum();
	}

    for(k = 0; k < MAX_SAMPLES && err < MAX_ERRORS; k++)
    {
		for(i = 1; i < tested; i++)
		{
			if (snd_out_d[k] != snd_out_d[i * MAX_SAMPLES + k])
			{
				printf("id:%d seg:%d sam:%d orig:%x here:%x\n", i, j, k, snd_out_d[k], snd_out_d[i * MAX_SAMPLES + k]);
				err++;
			}
		}
    }
  }

  printf("Test segments %d\n", j);
  for (k = 0; k < tested; k++)
  {
	  const char* myinfo = data_fn[k].fname;
	  printf("%d %s: %4.4f %4.4f %4.4f %4.4f\n", k, myinfo, elapsed[k], cputime[k], elapsed[k] / elapsed[0], cputime[k] / cputime[0]);
	  csv_put_float(elapsed[k]);
  }
  csv_put_string(",");
  for (k = 1; k < tested; k++)
  {
	  csv_put_float(elapsed[k] / elapsed[0]);
  }

  csv_put_string("\n");
  csv_close();

#endif
}
