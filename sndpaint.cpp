
#include <math.h>
#include <float.h>
#include <cassert>
#include <stdio.h>
#include <stdint.h>

#include "bdmpx.h"
#include "timing.h"
#include "platform.h"
#include "csv.h"
#include "intrin.h"

#define	MAX_QPATH		64

typedef float vec_t;
typedef vec_t vec2_t[2];
typedef vec_t vec3_t[3];
typedef vec_t vec4_t[4];
typedef vec_t vec5_t[5];

typedef enum { qfalse = 0, qtrue } qboolean;

#define	PAINTBUFFER_SIZE		4096					// this is in samples

#define SND_CHUNK_SIZE			1024					// samples
#define SND_CHUNK_SIZE_FLOAT	(SND_CHUNK_SIZE/2)		// floats
#define SND_CHUNK_SIZE_BYTE		(SND_CHUNK_SIZE*2)		// floats

typedef struct {
	int			left;	// the final values will be clamped to +/- 0x00ffff00 and shifted down
	int			right;
} portable_samplepair_t;

typedef struct adpcm_state {
	short	sample;		/* Previous output value */
	char	index;		/* Index into stepsize table */
} adpcm_state_t;

typedef	struct sndBuffer_s {
	short					sndChunk[SND_CHUNK_SIZE];
	struct sndBuffer_s		*next;
	int						size_notUsed;
	adpcm_state_t			adpcm_notUsed;
} sndBuffer;

typedef struct sfx_s {
	sndBuffer		*soundData;
	qboolean		defaultSound_notUsed;			// couldn't be loaded, so use buzz
	qboolean		inMemory_notUsed;				// not in Memory
	qboolean		soundCompressed_notUsed;		// not in Memory
	int				soundCompressionMethod_notUsed;
	int 			soundLength_notUsed;
	char 			soundName_notUsed[MAX_QPATH];
	int				lastTimeUsed_notUsed;
	struct sfx_s	*next_notUsed;
} sfx_t;

typedef struct
{
	int			allocTime_notUsed;
	int			startSample_notUsed;	// START_SAMPLE_IMMEDIATE = set immediately on next mix
	int			entnum_notUsed;			// to allow overriding a specific sound
	int			entchannel_notUsed;		// to allow overriding a specific sound
	int			leftvol;		// 0-255 volume after spatialization
	int			rightvol;		// 0-255 volume after spatialization
	int			master_vol_notUsed;		// 0-255 volume before spatialization
	float		dopplerScale;
	float		oldDopplerScale;
	vec3_t		origin_notUsed;			// only use if fixed_origin is set
	qboolean	fixed_origin_notUsed;	// use origin instead of fetching entnum's origin
	sfx_t		*thesfx_notUsed;		// sfx structure
	qboolean	doppler;
} channel_t;

static QALIGNA(16) portable_samplepair_t paintbuffer[2 * PAINTBUFFER_SIZE] QALIGNB(16);
static int snd_vol;

static void S_PaintChannelFrom16(channel_t *ch, const sfx_t *sc, int count, int sampleOffset, int bufferOffset) {
	int						data, aoff, boff;
	int						leftvol, rightvol;
	int						i, j;
	portable_samplepair_t	*samp;
	sndBuffer				*chunk;
	short					*samples;
	float					ooff, fdata, fdiv, fleftvol, frightvol;

	samp = &paintbuffer[bufferOffset];

	if (ch->doppler) {
		sampleOffset = sampleOffset * ch->oldDopplerScale;
	}

	chunk = sc->soundData;
	while (sampleOffset >= SND_CHUNK_SIZE) {
		chunk = chunk->next;
		sampleOffset -= SND_CHUNK_SIZE;
		if (!chunk) {
			chunk = sc->soundData;
		}
	}

	if (!ch->doppler || ch->dopplerScale == 1.0f) {
		leftvol = ch->leftvol*snd_vol;
		rightvol = ch->rightvol*snd_vol;
		samples = chunk->sndChunk;
		for (i = 0; i < count; i++) {
			data = samples[sampleOffset++];
			samp[i].left += (data * leftvol) >> 8;
			samp[i].right += (data * rightvol) >> 8;

			if (sampleOffset == SND_CHUNK_SIZE) {
				chunk = chunk->next;
				samples = chunk->sndChunk;
				sampleOffset = 0;
			}
		}
	}
	else {
		fleftvol = ch->leftvol*snd_vol;
		frightvol = ch->rightvol*snd_vol;

		ooff = sampleOffset;
		samples = chunk->sndChunk;

		for (i = 0; i < count; i++) {

			aoff = ooff;
			ooff = ooff + ch->dopplerScale;
			boff = ooff;
			fdata = 0;
			for (j = aoff; j < boff; j++) {
				if (j == SND_CHUNK_SIZE) {
					chunk = chunk->next;
					if (!chunk) {
						chunk = sc->soundData;
					}
					samples = chunk->sndChunk;
					ooff -= SND_CHUNK_SIZE;
				}
				fdata += samples[j&(SND_CHUNK_SIZE - 1)];
			}
			fdiv = 256 * (boff - aoff);
			samp[i].left += (fdata * fleftvol) / fdiv;
			samp[i].right += (fdata * frightvol) / fdiv;
		}
	}
}

#if idsse
#define SSEALIGNED(X) (((intptr_t)(X) & 15u) == 0)

//NOTE: _mm_mul_epi32 is SSE4.1, but SSE4.1 already has _mm_mullo_epi32 which does what we want
#define SSE_USE_MULLO 1
#if SSE_USE_MULLO
#define mymul_s32 _mm_mullo_epi32
#else
static ID_INLINE __m128i mymul_s32(const __m128i a, const __m128i b)
{
	__m128i tmp1 = _mm_mul_epi32(a, b); /* mul 2,0*/
	__m128i tmp2 = _mm_mul_epi32(_mm_srli_si128(a, 4), _mm_srli_si128(b, 4)); /* mul 3,1 */
	return _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(tmp1), _mm_castsi128_ps(tmp2), _MM_SHUFFLE(2, 0, 2, 0))); /* shuffle results and pack */
}
#endif
static void S_PaintChannelFrom16_SSE(channel_t *ch, const sfx_t *sc, int count, int sampleOffset, int bufferOffset) {
	int						data, aoff, boff;
	int						leftvol, rightvol;
	int						i, j;
	portable_samplepair_t	*samp;
	sndBuffer				*chunk;
	short					*samples;
	float					ooff, fdata, fdiv, fleftvol, frightvol;

	samp = &paintbuffer[bufferOffset];

	if (ch->doppler) {
		sampleOffset = sampleOffset * ch->oldDopplerScale;
	}

	chunk = sc->soundData;
	while (sampleOffset >= SND_CHUNK_SIZE) {
		chunk = chunk->next;
		sampleOffset -= SND_CHUNK_SIZE;
		if (!chunk) {
			chunk = sc->soundData;
		}
	}

	if (!ch->doppler || ch->dopplerScale == 1.0f) {
		__m128i volume_vec;
		int vectorCount, samplesLeft, chunkSamplesLeft;

		__m128i s0_s16;
		__m128i sampleData_a, sampleData_b;
		__m128i sampleData0, sampleData1, sampleData2, sampleData3;
		__m128i d0_s32, d1_s32, d2_s32, d3_s32;

		leftvol = ch->leftvol*snd_vol;
		rightvol = ch->rightvol*snd_vol;
		samples = chunk->sndChunk;

#if SSE_USE_MULLO
		volume_vec = _mm_setr_epi32(leftvol, rightvol, leftvol, rightvol);
#else
		volume_vec = _mm_setr_epi32(leftvol, leftvol, rightvol, rightvol);
#endif

		i = 0;

		while (i < count) {
			/* Try to align destination to 16-byte boundary */
			while (i < count && (!SSEALIGNED(&samp[i]) || ((count - i) < 8) || ((SND_CHUNK_SIZE - sampleOffset) < 8))) {
				data = samples[sampleOffset++];
				samp[i].left += (data * leftvol) >> 8;
				samp[i].right += (data * rightvol) >> 8;

				if (sampleOffset == SND_CHUNK_SIZE) {
					chunk = chunk->next;
					samples = chunk->sndChunk;
					sampleOffset = 0;
				}
				i++;
			}
			/* Destination is now aligned.  Process as many 8-sample
			chunks as we can before we run out of room from the current
			sound chunk.  We do 8 per loop to avoid extra source data reads. */
			samplesLeft = count - i;
			chunkSamplesLeft = SND_CHUNK_SIZE - sampleOffset;
			if (samplesLeft > chunkSamplesLeft)
				samplesLeft = chunkSamplesLeft;

			vectorCount = samplesLeft / 8;

			while (vectorCount)
			{
				//s0_s16 = _mm_loadu_si128((__m128i const*)&samples[sampleOffset]);
				s0_s16 = _mm_lddqu_si128((__m128i const*)&samples[sampleOffset]); //SSE3

				sampleData_a = _mm_unpacklo_epi16(_mm_setzero_si128(), s0_s16);
				sampleData_b = _mm_unpackhi_epi16(_mm_setzero_si128(), s0_s16);

#if SSE_USE_MULLO
				sampleData0 = _mm_shuffle_epi32(sampleData_a, _MM_SHUFFLE(1, 1, 0, 0));
				sampleData1 = _mm_shuffle_epi32(sampleData_a, _MM_SHUFFLE(3, 3, 2, 2));
#else
				sampleData0 = _mm_shuffle_epi32(sampleData_a, _MM_SHUFFLE(1, 0, 1, 0));
				sampleData1 = _mm_shuffle_epi32(sampleData_a, _MM_SHUFFLE(3, 2, 3, 2));
#endif

				sampleData0 = _mm_srai_epi32(sampleData0, 16);
				sampleData0 = mymul_s32(sampleData0, volume_vec);
				sampleData0 = _mm_srai_epi32(sampleData0, 8);

				sampleData1 = _mm_srai_epi32(sampleData1, 16);
				sampleData1 = mymul_s32(sampleData1, volume_vec);
				sampleData1 = _mm_srai_epi32(sampleData1, 8);

				/* Load up destination sample data */
				d0_s32 = _mm_load_si128((__m128i*)&samp[i]);
				d1_s32 = _mm_load_si128((__m128i*)&samp[i + 2]);
				d2_s32 = _mm_load_si128((__m128i*)&samp[i + 4]);
				d3_s32 = _mm_load_si128((__m128i*)&samp[i + 6]);

				d0_s32 = _mm_add_epi32(d0_s32, sampleData0);
				d1_s32 = _mm_add_epi32(d1_s32, sampleData1);


#if SSE_USE_MULLO
				sampleData2 = _mm_shuffle_epi32(sampleData_b, _MM_SHUFFLE(1, 1, 0, 0));
				sampleData3 = _mm_shuffle_epi32(sampleData_b, _MM_SHUFFLE(3, 3, 2, 2));
#else
				sampleData2 = _mm_shuffle_epi32(sampleData_b, _MM_SHUFFLE(1, 0, 1, 0));
				sampleData3 = _mm_shuffle_epi32(sampleData_b, _MM_SHUFFLE(3, 2, 3, 2));
#endif

				sampleData2 = _mm_srai_epi32(sampleData2, 16);
				sampleData2 = mymul_s32(sampleData2, volume_vec);
				sampleData2 = _mm_srai_epi32(sampleData2, 8);

				sampleData3 = _mm_srai_epi32(sampleData3, 16);
				sampleData3 = mymul_s32(sampleData3, volume_vec);
				sampleData3 = _mm_srai_epi32(sampleData3, 8);

				d2_s32 = _mm_add_epi32(d2_s32, sampleData2);
				d3_s32 = _mm_add_epi32(d3_s32, sampleData3);

				/* Store destination sample data */
				_mm_store_si128((__m128i*)&samp[i], d0_s32);
				_mm_store_si128((__m128i*)&samp[i + 2], d1_s32);
				_mm_store_si128((__m128i*)&samp[i + 4], d2_s32);
				_mm_store_si128((__m128i*)&samp[i + 6], d3_s32);

				i += 8;
				vectorCount--;
				sampleOffset += 8;

				if (sampleOffset == SND_CHUNK_SIZE) {
					chunk = chunk->next;
					samples = chunk->sndChunk;
					sampleOffset = 0;
				}
			}
		}
	}
	else {
		fleftvol = ch->leftvol*snd_vol;
		frightvol = ch->rightvol*snd_vol;

		ooff = sampleOffset;
		samples = chunk->sndChunk;

		for (i = 0; i < count; i++) {

			aoff = ooff;
			ooff = ooff + ch->dopplerScale;
			boff = ooff;
			fdata = 0;
			for (j = aoff; j < boff; j++) {
				if (j == SND_CHUNK_SIZE) {
					chunk = chunk->next;
					if (!chunk) {
						chunk = sc->soundData;
					}
					samples = chunk->sndChunk;
					ooff -= SND_CHUNK_SIZE;
				}
				fdata += samples[j&(SND_CHUNK_SIZE - 1)];
			}
			fdiv = 256 * (boff - aoff);
			samp[i].left += (fdata * fleftvol) / fdiv;
			samp[i].right += (fdata * frightvol) / fdiv;
		}
	}
}
#endif

#if idneon
static void S_PaintChannelFrom16_NEON(channel_t *ch, const sfx_t *sc, int count, int sampleOffset, int bufferOffset) {
	int						data, aoff, boff;
	int						leftvol, rightvol;
	int						i, j;
	portable_samplepair_t	*samp;
	sndBuffer				*chunk;
	short					*samples;
	float					ooff, fdata, fdiv, fleftvol, frightvol;

	samp = &paintbuffer[bufferOffset];

	if (ch->doppler) {
		sampleOffset = sampleOffset * ch->oldDopplerScale;
	}

	chunk = sc->soundData;
	while (sampleOffset >= SND_CHUNK_SIZE) {
		chunk = chunk->next;
		sampleOffset -= SND_CHUNK_SIZE;
		if (!chunk) {
			chunk = sc->soundData;
		}
	}

	if (!ch->doppler || ch->dopplerScale == 1.0f) {
		int32x4_t volume_vec, snd_vol_vec;
		int vectorCount, samplesLeft, chunkSamplesLeft;

		snd_vol_vec = vld1q_dup_s32(&snd_vol);
		int32x2_t volume_vec_d0 = vld1_s32((int32_t*)&ch->leftvol); //load leftvol and rightvol
		int32x2_t volume_vec_d1 = volume_vec_d0;

		volume_vec = vcombine_s32(volume_vec_d0, volume_vec_d1);
		volume_vec = vmulq_s32(volume_vec, snd_vol_vec);

		leftvol = ch->leftvol*snd_vol;
		rightvol = ch->rightvol*snd_vol;
		samples = chunk->sndChunk;

		i = 0;

		while (i < count) {
			/* Try to align destination to 16-byte boundary */
			while (i < count && (((unsigned long)&samp[i] & 15u) || ((count - i) < 8) || ((SND_CHUNK_SIZE - sampleOffset) < 8))) {
				data = samples[sampleOffset++];
				samp[i].left += (data * leftvol) >> 8;
				samp[i].right += (data * rightvol) >> 8;

				if (sampleOffset == SND_CHUNK_SIZE) {
					chunk = chunk->next;
					samples = chunk->sndChunk;
					sampleOffset = 0;
				}
				i++;
			}
			/* Destination is now aligned.  Process as many 8-sample
			chunks as we can before we run out of room from the current
			sound chunk.  We do 8 per loop to avoid extra source data reads. */
			samplesLeft = count - i;
			chunkSamplesLeft = SND_CHUNK_SIZE - sampleOffset;
			if (samplesLeft > chunkSamplesLeft)
				samplesLeft = chunkSamplesLeft;

			vectorCount = samplesLeft / 8;

			while (vectorCount)
			{
				int16x8_t samples_short = vld1q_s16(&samples[sampleOffset]);
				uint32x4_t samples_int0 = (uint32x4_t)vmovl_s16(vget_low_s16(samples_short));

				uint64x2_t samples_0 = vmovl_u32(vget_low_u32(samples_int0));
				uint64x2_t samples_0t = vshlq_n_u64(samples_0, 32);
				samples_0 = vorrq_u64(samples_0, samples_0t);

				uint64x2_t samples_1 = vmovl_u32(vget_high_u32(samples_int0));
				uint64x2_t samples_1t = vshlq_n_u64(samples_1, 32);
				samples_1 = vorrq_u64(samples_1, samples_1t);

				uint32x4_t samples_int1 = (uint32x4_t)vmovl_s16(vget_high_s16(samples_short));

				uint64x2_t samples_2 = vmovl_u32(vget_low_u32(samples_int1));
				uint64x2_t samples_2t = vshlq_n_u64(samples_2, 32);
				samples_2 = vorrq_u64(samples_2, samples_2t);

				uint64x2_t samples_3 = vmovl_u32(vget_high_u32(samples_int1));
				uint64x2_t samples_3t = vshlq_n_u64(samples_3, 32);
				samples_3 = vorrq_u64(samples_3, samples_3t);

				int32x4_t d0_s32 = vld1q_s32((int32_t*)&samp[i]);
				int32x4_t d1_s32 = vld1q_s32((int32_t*)&samp[i + 2]);
				int32x4_t d2_s32 = vld1q_s32((int32_t*)&samp[i + 3]);
				int32x4_t d3_s32 = vld1q_s32((int32_t*)&samp[i + 4]);

				int32x4_t merge0_s32 = vmulq_s32((int32x4_t)samples_0, volume_vec);
				merge0_s32 = vshrq_n_s32(merge0_s32, 8);

				int32x4_t merge1_s32 = vmulq_s32((int32x4_t)samples_1, volume_vec);
				merge1_s32 = vshrq_n_s32(merge1_s32, 8);

				d0_s32 = vaddq_s32(d0_s32, merge0_s32);
				d1_s32 = vaddq_s32(d1_s32, merge1_s32);


				int32x4_t merge2_s32 = vmulq_s32((int32x4_t)samples_2, volume_vec);
				merge2_s32 = vshrq_n_s32(merge2_s32, 8);

				int32x4_t merge3_s32 = vmulq_s32((int32x4_t)samples_3, volume_vec);
				merge3_s32 = vshrq_n_s32(merge3_s32, 8);

				d2_s32 = vaddq_s32(d2_s32, merge2_s32);
				d3_s32 = vaddq_s32(d3_s32, merge3_s32);

				vst1q_s32((int32_t*)&samp[i], d0_s32);
				vst1q_s32((int32_t*)&samp[i + 2], d1_s32);
				vst1q_s32((int32_t*)&samp[i + 4], d2_s32);
				vst1q_s32((int32_t*)&samp[i + 6], d3_s32);

				i += 8;
				vectorCount--;
				sampleOffset += 8;

				if (sampleOffset == SND_CHUNK_SIZE) {
					chunk = chunk->next;
					samples = chunk->sndChunk;
					sampleOffset = 0;
				}
			}
		}
	}
	else {
		fleftvol = ch->leftvol*snd_vol;
		frightvol = ch->rightvol*snd_vol;

		ooff = sampleOffset;
		samples = chunk->sndChunk;

		for (i = 0; i < count; i++) {

			aoff = ooff;
			ooff = ooff + ch->dopplerScale;
			boff = ooff;
			fdata = 0;
			for (j = aoff; j < boff; j++) {
				if (j == SND_CHUNK_SIZE) {
					chunk = chunk->next;
					if (!chunk) {
						chunk = sc->soundData;
					}
					samples = chunk->sndChunk;
					ooff -= SND_CHUNK_SIZE;
				}
				fdata += samples[j&(SND_CHUNK_SIZE - 1)];
			}
			fdiv = 256 * (boff - aoff);
			samp[i].left += (fdata * fleftvol) / fdiv;
			samp[i].right += (fdata * frightvol) / fdiv;
		}
	}
}
#endif

#if idarm
static void S_PaintChannelFrom16_ARMv6(channel_t *ch, const sfx_t *sc, int count, int sampleOffset, int bufferOffset) {
	int						aoff, boff;
	unsigned int					datau;
	int						leftvol, rightvol;
	int						i, j;
	portable_samplepair_t	*samp;
	sndBuffer				*chunk;
	short						*samples;
	float					ooff, fdata, fdiv, fleftvol, frightvol;

	samp = &paintbuffer[bufferOffset];

	if (ch->doppler) {
		sampleOffset = sampleOffset * ch->oldDopplerScale;
	}

	chunk = sc->soundData;
	while (sampleOffset >= SND_CHUNK_SIZE) {
		chunk = chunk->next;
		sampleOffset -= SND_CHUNK_SIZE;
		if (!chunk) {
			chunk = sc->soundData;
		}
	}

	if (!ch->doppler || ch->dopplerScale == 1.0f) {
		int samp_l, samp_r, samp_l1, samp_r1;
		leftvol = (ch->leftvol*snd_vol) << 8; //will get shifted >> 16 later on instead of >> 8, so we compensate here
		rightvol = (ch->rightvol*snd_vol) << 8;
		samples = chunk->sndChunk;
		for (i = 0; i < count; i += 2) {
			//NOTE: should we check for samples being aligned on 4 bytes?
			datau = *(unsigned int*)(&samples[sampleOffset]);
			sampleOffset += 2;

			samp_l = samp[i].left;
			samp_r = samp[i].right;
			asm("smlawb %0, %1, %2, %3" : "+r" (samp_l) : "r" (leftvol), "r" (datau), "r" (samp_l));
			asm("smlawb %0, %1, %2, %3" : "+r" (samp_r) : "r" (rightvol), "r" (datau), "r" (samp_r));
			samp[i].left = samp_l;
			samp[i].right = samp_r;

			samp_l1 = samp[i + 1].left;
			samp_r1 = samp[i + 1].right;
			asm("smlawt %0, %1, %2, %3" : "+r" (samp_l1) : "r" (leftvol), "r" (datau), "r" (samp_l1));
			asm("smlawt %0, %1, %2, %3" : "+r" (samp_r1) : "r" (rightvol), "r" (datau), "r" (samp_r1));
			samp[i + 1].left = samp_l1;
			samp[i + 1].right = samp_r1;

			if (sampleOffset == SND_CHUNK_SIZE) {
				chunk = chunk->next;
				samples = chunk->sndChunk;
				sampleOffset = 0;
			}
		}
	}
	else {
		fleftvol = ch->leftvol*snd_vol;
		frightvol = ch->rightvol*snd_vol;

		ooff = sampleOffset;
		samples = chunk->sndChunk;

		for (i = 0; i < count; i++) {

			aoff = ooff;
			ooff = ooff + ch->dopplerScale;
			boff = ooff;
			fdata = 0;
			for (j = aoff; j < boff; j++) {
				if (j == SND_CHUNK_SIZE) {
					chunk = chunk->next;
					if (!chunk) {
						chunk = sc->soundData;
					}
					samples = chunk->sndChunk;
					ooff -= SND_CHUNK_SIZE;
				}
				fdata += samples[j&(SND_CHUNK_SIZE - 1)];
			}
			fdiv = 256 * (boff - aoff);
			samp[i].left += (fdata * fleftvol) / fdiv;
			samp[i].right += (fdata * frightvol) / fdiv;
		}
	}
}
#endif

#if idneon
void myfunc(channel_t *ch, short *samples)
{
	int i = 0;
	portable_samplepair_t	*samp;

	samp = &paintbuffer[0];

	int32x4_t volume_vec, snd_vol_vec;
	int vectorCount, samplesLeft, chunkSamplesLeft;

	snd_vol_vec = vld1q_dup_s32(&snd_vol);
	int32x2_t volume_vec_d0 = vld1_s32((int32_t*)&ch->leftvol); //load leftvol and rightvol
	int32x2_t volume_vec_d1 = volume_vec_d0;

	volume_vec = vcombine_s32(volume_vec_d0, volume_vec_d1);
	volume_vec = vmulq_s32(volume_vec, snd_vol_vec);


	int16x8_t samples_short = vld1q_s16(&samples[0]);
	uint32x4_t samples_int0 = (uint32x4_t)vmovl_s16(vget_low_s16(samples_short));

	uint64x2_t samples_0 = vmovl_u32(vget_low_u32(samples_int0));
	uint64x2_t samples_0t = vshlq_n_u64(samples_0, 32);
	samples_0 = vorrq_u64(samples_0, samples_0t);

	uint64x2_t samples_1 = vmovl_u32(vget_high_u32(samples_int0));
	uint64x2_t samples_1t = vshlq_n_u64(samples_1, 32);
	samples_1 = vorrq_u64(samples_1, samples_1t);

	uint32x4_t samples_int1 = (uint32x4_t)vmovl_s16(vget_high_s16(samples_short));

	uint64x2_t samples_2 = vmovl_u32(vget_low_u32(samples_int1));
	uint64x2_t samples_2t = vshlq_n_u64(samples_2, 32);
	samples_2 = vorrq_u64(samples_2, samples_2t);

	uint64x2_t samples_3 = vmovl_u32(vget_high_u32(samples_int1));
	uint64x2_t samples_3t = vshlq_n_u64(samples_3, 32);
	samples_3 = vorrq_u64(samples_3, samples_3t);

	int32x4_t merge0_s32 = vmulq_s32((int32x4_t)samples_0, volume_vec);
	merge0_s32 = vshrq_n_s32(merge0_s32, 8);

	int32x4_t merge1_s32 = vmulq_s32((int32x4_t)samples_1, volume_vec);
	merge1_s32 = vshrq_n_s32(merge1_s32, 8);


	int32x4_t merge2_s32 = vmulq_s32((int32x4_t)samples_2, volume_vec);
	merge2_s32 = vshrq_n_s32(merge2_s32, 8);

	int32x4_t merge3_s32 = vmulq_s32((int32x4_t)samples_3, volume_vec);
	merge3_s32 = vshrq_n_s32(merge3_s32, 8);


	int32x4_t d0_s32 = vld1q_s32((int32_t*)&samp[i]);
	int32x4_t d1_s32 = vld1q_s32((int32_t*)&samp[i + 2]);
	int32x4_t d2_s32 = vld1q_s32((int32_t*)&samp[i + 3]);
	int32x4_t d3_s32 = vld1q_s32((int32_t*)&samp[i + 4]);

	d0_s32 = vaddq_s32(d0_s32, merge0_s32);
	d1_s32 = vaddq_s32(d1_s32, merge1_s32);

	d2_s32 = vaddq_s32(d2_s32, merge2_s32);
	d3_s32 = vaddq_s32(d3_s32, merge3_s32);

	vst1q_s32((int32_t*)&samp[i], d0_s32);
	vst1q_s32((int32_t*)&samp[i + 2], d1_s32);
	vst1q_s32((int32_t*)&samp[i + 4], d2_s32);
	vst1q_s32((int32_t*)&samp[i + 6], d3_s32);
}
#endif

typedef void(*sndpaint_fn)(channel_t *ch, const sfx_t *sc, int count, int sampleOffset, int bufferOffset);

struct function_data
{
	void* fp;
	const char* fname;
};

static struct function_data data_fn[] =
{
	{ S_PaintChannelFrom16,        "      sndpaint" },
#if idsse
	{ S_PaintChannelFrom16_SSE,    "  sndpaint_sse" },
#endif
#if idneon
	{ S_PaintChannelFrom16_NEON,   " sndpaint_neon" },
#endif
#if idarm
	{ S_PaintChannelFrom16_ARMv6,  "sndpaint_armv6" },
#endif
	{ 0, 0 }
};

static sndBuffer mychunks[10];

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(X[0]))

void myfunc(channel_t *ch, short *samples);

#define REPETITIONS 1
#define MAX_ERRORS 10
#define TEST_SAMPLESIZE 900000
void maintest_sndpaint(void)
{
#if tests
	channel_t ch;

	snd_vol = 4253;

	ch.leftvol = 1;
	ch.rightvol = 1;

	short samples[8] = { -1, -2, -3, -4, -5, -6, -7, -8 };

	myfunc(&ch, samples);

	printf("%8d %8d %8d %8d\n", paintbuffer[0].left, paintbuffer[0].right, paintbuffer[1].left, paintbuffer[1].right);
	printf("%8d %8d %8d %8d\n", paintbuffer[2].left, paintbuffer[2].right, paintbuffer[3].left, paintbuffer[3].right);
	printf("%8d %8d %8d %8d\n", paintbuffer[4].left, paintbuffer[4].right, paintbuffer[5].left, paintbuffer[5].right);
	printf("%8d %8d %8d %8d\n", paintbuffer[6].left, paintbuffer[6].right, paintbuffer[7].left, paintbuffer[7].right);
#endif

	int i, j, k;
	int err = 0;
	int ercd;
	int tested = 0;
	Timer timer[ARRAY_SIZE(data_fn) - 1];
	Cputime cpu[ARRAY_SIZE(data_fn) - 1];
	double elapsed[ARRAY_SIZE(data_fn) - 1];
	double cputime[ARRAY_SIZE(data_fn) - 1];

	Timer progtime;

	if (ARRAY_SIZE(data_fn) - 1 > ARRAY_SIZE(paintbuffer) / PAINTBUFFER_SIZE)
	{
		printf("Need to increase size of paintbuffer to: %d * PAINTBUFFER_SIZE\n", (int)ARRAY_SIZE(data_fn) - 1);
		return;
	}

	ercd = csv_open("./results_sndpaint.csv");
	if (ercd != 0)
	{
		printf("Could not open the csv file.\n\n");
	}

	//bdmpx_set_option(BDMPX_OPTION_PRECACHE_FILE_FOR_READ, 1);
	ercd = bdmpx_create(NULL, "bdmpx_paintchannel.bin", BDMPX_OP_READ);
	if (ercd != 0)
	{
		printf("Could not open the test input file. Aborting.\n");
		return;
	}

	int32_t count;
	int32_t sampleOffset;
	channel_t channel;
	sfx_t sfx;
	const char assert_int_size[sizeof(int32_t) == sizeof(int)] = { 0 };
	int total_processed = 0;

	snd_vol = 204;

	memset(&sfx, 0, sizeof(sfx));
	sfx.soundData = mychunks;
	for (i = 0; i < ARRAY_SIZE(mychunks); i++)
	{
		mychunks[i].next = &mychunks[i + 1];
	}

#if 0
	for (int x = 0; x < 100; x++)
	{
		int a, b, c, d, e;
		a = b = c = d = e = INT_MAX;
		ercd = bdmpx_read(NULL, 5, &a, 0, &b, 0, &c, 0, &d, 0, &e, 0);
		printf("%d: %d %d %d %d %d\n", ercd, a, b, c, d, e);
	}
	return;
#endif

	for (j = 0; j < TEST_SAMPLESIZE && err < MAX_ERRORS; j++)
	{
		int chunkidx = 0;
		ercd = bdmpx_read(NULL, 3, 0, &count, 0, &channel, 0, &sampleOffset, 0, &snd_vol);
		if (ercd != 3)
		{
			break;
		}
		while ((SND_CHUNK_SIZE * chunkidx < sampleOffset + count) || (1 == bdmpx_peek(NULL)))
		{
			if (chunkidx > ARRAY_SIZE(mychunks)) { printf("too many chunks\n"); return; }
			int chnksz = sizeof(mychunks[0].sndChunk);
			ercd = bdmpx_read(NULL, 1, &chnksz, &mychunks[chunkidx]);
			if (chnksz < sizeof(mychunks[0].sndChunk)) { printf("less than expected\n"); };
			if (ercd != 1) { printf("reading chunks error, found:%d elems\n", ercd); return; }
			chunkidx++;
		}
		if (sampleOffset + count > SND_CHUNK_SIZE * chunkidx)
		{
			printf("too few chunks: off:%d count:%d read:%d needed:%d\n", sampleOffset, count, SND_CHUNK_SIZE * chunkidx, sampleOffset + count);
			return;
		}

		total_processed += count;
		memset(paintbuffer, 0, sizeof(paintbuffer));

		for (tested = 0; data_fn[tested].fp != 0; tested++)
		{
			sndpaint_fn callme = (sndpaint_fn)data_fn[tested].fp;
			const char* myinfo = data_fn[tested].fname;
			
			cpu[tested].reset();
			timer[tested].reset();
			for (int r = 0; r < REPETITIONS; r++)
				callme(&channel, &sfx, count, sampleOffset, tested * PAINTBUFFER_SIZE);

			elapsed[tested] = timer[tested].accum();
			cputime[tested] = cpu[tested].accum();
		}

		for (i = 1; i < tested; i++)
		{
			for (k = 0; k < count && err < MAX_ERRORS; k++)
			{
				portable_samplepair_t	*samp0 = &paintbuffer[0];
				portable_samplepair_t	*samp1 = &paintbuffer[i * PAINTBUFFER_SIZE];
				if (samp0[k].left != samp1[k].left || samp0[k].right != samp1[k].right)
				{
					printf("id:%d seg:%d sam:%d origL:%8x hereL:%8x  origR:%8x hereR:%8x\n", i, j, k, samp0[k].left, samp1[k].left, samp0[k].right, samp1[k].right);
					err++;
				}
			}
		}

	}

	printf("SndPaint Test segments %d, total samples %d\n", j, total_processed);
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

	printf("SndPaint Prog time: %4.4f sec\n", progtime.elapsed_ms() / 1000);

	printf("\n");
}