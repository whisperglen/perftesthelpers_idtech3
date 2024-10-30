
#include <math.h>
#include <float.h>
#include <cassert>
#include <intrin.h>
#include <stdio.h>

#include "bdmpx.h"
#include "timing.h"
#include "platform.h"

#define	myftol(x) ((int)(x))
#define Q_fabs fabsf

#define DotProduct(x,y)			((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2])
#define VectorSubtract(a,b,c)	((c)[0]=(a)[0]-(b)[0],(c)[1]=(a)[1]-(b)[1],(c)[2]=(a)[2]-(b)[2])
#define VectorAdd(a,b,c)		((c)[0]=(a)[0]+(b)[0],(c)[1]=(a)[1]+(b)[1],(c)[2]=(a)[2]+(b)[2])
#define VectorCopy(a,b)			((b)[0]=(a)[0],(b)[1]=(a)[1],(b)[2]=(a)[2])
#define	VectorScale(v, s, o)	((o)[0]=(v)[0]*(s),(o)[1]=(v)[1]*(s),(o)[2]=(v)[2]*(s))
#define	VectorMA(v, s, b, o)	((o)[0]=(v)[0]+(b)[0]*(s),(o)[1]=(v)[1]+(b)[1]*(s),(o)[2]=(v)[2]+(b)[2]*(s))

#define LERP( a, b, w ) ( ( a ) * ( 1.0f - ( w ) ) + ( b ) * ( w ) )
#define LUMA( red, green, blue ) ( 0.2126f * ( red ) + 0.7152f * ( green ) + 0.0722f * ( blue ) )

typedef int		qhandle_t;
typedef unsigned short glIndex_t;

typedef float vec_t;
typedef vec_t vec2_t[2];
typedef vec_t vec3_t[3];
typedef vec_t vec4_t[4];

typedef enum { qfalse = 0, qtrue } qboolean;

// trRefdef_t holds everything that comes in refdef_t,
// as well as the locally generated scene information
typedef struct {
	int			num_dlights;
	struct dlight_s	*dlights;
} trRefdef_t;

typedef struct {
	int		c_surfaces, c_shaders, c_vertexes, c_indexes, c_totalIndexes;
	float	c_overDraw;

	int		c_dlightVertexes;
	int		c_dlightIndexes;

	int		c_flareAdds;
	int		c_flareTests;
	int		c_flareRenders;

	int		msec;			// total msec for backend run
} backEndCounters_t;

static struct {
	trRefdef_t	refdef;
	backEndCounters_t	pc;
} backEnd;

typedef struct dlight_s {
	vec3_t	origin;
	vec3_t	color;				// range from 0.0 to 1.0, should be color normalized
	float	radius;

	vec3_t	transformed;		// origin in local coordinate system
	int		additive;			// texture detail is lost tho when the lightmap is dark
} dlight_t;


// surface geometry should not exceed these limits
#define	SHADER_MAX_VERTEXES	1000
#define	SHADER_MAX_INDEXES	(6*SHADER_MAX_VERTEXES)

static struct
{
	int			dlightBits;	// or together of all vertexDlightBits
	int         numVertexes;
	int			numIndexes;
	QALIGNA(16) glIndex_t	indexes[SHADER_MAX_INDEXES];
	QALIGNA(16) vec4_t		xyz[SHADER_MAX_VERTEXES];
	QALIGNA(16) vec4_t		normal[SHADER_MAX_VERTEXES];
} tess;

static struct rval
{
	int integer;
	float value;
} greyscale = { 0,0. }, dlightBacks = { 1,1. }, *r_greyscale = &greyscale, *r_dlightBacks = &dlightBacks;

static void ProjectDlightTexture_scalar(void) {
	int		i, l;
	vec3_t	origin;
	float	*texCoords;
	byte	*colors;
	byte	clipBits[SHADER_MAX_VERTEXES];
	float	texCoordsArray[SHADER_MAX_VERTEXES][2];
	byte	colorArray[SHADER_MAX_VERTEXES][4];
	glIndex_t	hitIndexes[SHADER_MAX_INDEXES];
	int		numIndexes;
	float	scale;
	float	radius;
	vec3_t	floatColor;
	float	modulate = 0.0f;

	if (!backEnd.refdef.num_dlights) {
		return;
	}

	for (l = 0; l < backEnd.refdef.num_dlights; l++) {
		dlight_t	*dl;

		if (!(tess.dlightBits & (1 << l))) {
			continue;	// this surface definately doesn't have any of this light
		}
		texCoords = texCoordsArray[0];
		colors = colorArray[0];

		dl = &backEnd.refdef.dlights[l];
		VectorCopy(dl->transformed, origin);
		radius = dl->radius;
		scale = 1.0f / radius;

		if (r_greyscale->integer)
		{
			float luminance;

			luminance = LUMA(dl->color[0], dl->color[1], dl->color[2]) * 255.0f;
			floatColor[0] = floatColor[1] = floatColor[2] = luminance;
		}
		else if (r_greyscale->value)
		{
			float luminance;

			luminance = LUMA(dl->color[0], dl->color[1], dl->color[2]) * 255.0f;
			floatColor[0] = LERP(dl->color[0] * 255.0f, luminance, r_greyscale->value);
			floatColor[1] = LERP(dl->color[1] * 255.0f, luminance, r_greyscale->value);
			floatColor[2] = LERP(dl->color[2] * 255.0f, luminance, r_greyscale->value);
		}
		else
		{
			floatColor[0] = dl->color[0] * 255.0f;
			floatColor[1] = dl->color[1] * 255.0f;
			floatColor[2] = dl->color[2] * 255.0f;
		}

		for (i = 0; i < tess.numVertexes; i++, texCoords += 2, colors += 4) {
			int		clip = 0;
			vec3_t	dist;

			VectorSubtract(origin, tess.xyz[i], dist);

			backEnd.pc.c_dlightVertexes++;

			texCoords[0] = 0.5f + dist[0] * scale;
			texCoords[1] = 0.5f + dist[1] * scale;

			if (!r_dlightBacks->integer &&
				// dist . tess.normal[i]
				(dist[0] * tess.normal[i][0] +
					dist[1] * tess.normal[i][1] +
					dist[2] * tess.normal[i][2]) < 0.0f) {
				clip = 63;
			}
			else {
				if (texCoords[0] < 0.0f) {
					clip |= 1;
				}
				else if (texCoords[0] > 1.0f) {
					clip |= 2;
				}
				if (texCoords[1] < 0.0f) {
					clip |= 4;
				}
				else if (texCoords[1] > 1.0f) {
					clip |= 8;
				}

				// modulate the strength based on the height and color
				if (dist[2] > radius) {
					clip |= 16;
					modulate = 0.0f;
				}
				else if (dist[2] < -radius) {
					clip |= 32;
					modulate = 0.0f;
				}
				else {
					dist[2] = Q_fabs(dist[2]);
					if (dist[2] < radius * 0.5f) {
						modulate = 1.0f;
					}
					else {
						modulate = 2.0f * (radius - dist[2]) * scale;
					}
				}
			}
			clipBits[i] = clip;
			colors[0] = myftol(floatColor[0] * modulate);
			colors[1] = myftol(floatColor[1] * modulate);
			colors[2] = myftol(floatColor[2] * modulate);
			colors[3] = 255;
		}

		// build a list of triangles that need light
		numIndexes = 0;
		for (i = 0; i < tess.numIndexes; i += 3) {
			glIndex_t	a, b, c;

			a = tess.indexes[i];
			b = tess.indexes[i + 1];
			c = tess.indexes[i + 2];
			if (clipBits[a] & clipBits[b] & clipBits[c]) {
				continue;	// not lighted
			}
			hitIndexes[numIndexes] = a;
			hitIndexes[numIndexes + 1] = b;
			hitIndexes[numIndexes + 2] = c;
			numIndexes += 3;
		}

		if (!numIndexes) {
			continue;
		}
#if 0
		qglEnableClientState(GL_TEXTURE_COORD_ARRAY);
		qglTexCoordPointer(2, GL_FLOAT, 0, texCoordsArray[0]);

		qglEnableClientState(GL_COLOR_ARRAY);
		qglColorPointer(4, GL_UNSIGNED_BYTE, 0, colorArray);

		GL_Bind(tr.dlightImage);
		// include GLS_DEPTHFUNC_EQUAL so alpha tested surfaces don't add light
		// where they aren't rendered
		if (dl->additive) {
			GL_State(GLS_SRCBLEND_ONE | GLS_DSTBLEND_ONE | GLS_DEPTHFUNC_EQUAL);
		}
		else {
			GL_State(GLS_SRCBLEND_DST_COLOR | GLS_DSTBLEND_ONE | GLS_DEPTHFUNC_EQUAL);
		}
		R_DrawElements(numIndexes, hitIndexes);
#endif
		backEnd.pc.c_totalIndexes += numIndexes;
		backEnd.pc.c_dlightIndexes += numIndexes;
	}
}

static void ProjectDlightTexture_vector(void) {
	int		i, l;
	__m128 modulateVec, floatColorVec;
	__m128 vSel, vSel1;
	__m128i shuf;
	vec3_t	origin;
	float	*texCoords;
	byte	*colors;
	byte	clipBits[SHADER_MAX_VERTEXES];
	float	texCoordsArray[SHADER_MAX_VERTEXES][2];
	QALIGNA(16) byte	colorArray[SHADER_MAX_VERTEXES][4];
	glIndex_t	hitIndexes[SHADER_MAX_INDEXES];
	int		numIndexes;
	float	scale;
	float	radius;
	QALIGNA(16) vec4_t	floatColor;
	float	modulate = 0.0f;

	if (!backEnd.refdef.num_dlights) {
		return;
	}

	floatColor[3] = 1.0f;

	vSel = _mm_set1_ps(255);
	vSel1 = _mm_setr_ps(0, 0, 0, 255);

	shuf = _mm_setr_epi8(3, 7, 11, 15, 2, 6, 10, 14, 1, 5, 9, 13, 0, 4, 8, 12);

	for (l = 0; l < backEnd.refdef.num_dlights; l++) {
		dlight_t	*dl;
		__m128i colorVec = _mm_setzero_si128();
		int step = 0;

		if (!(tess.dlightBits & (1 << l))) {
			continue;	// this surface definately doesn't have any of this light
		}
		texCoords = texCoordsArray[0];
		colors = colorArray[0];

		dl = &backEnd.refdef.dlights[l];
		VectorCopy(dl->transformed, origin);
		radius = dl->radius;
		scale = 1.0f / radius;

		if (r_greyscale->integer)
		{
			float luminance;

			luminance = LUMA(dl->color[0], dl->color[1], dl->color[2]) * 255.0f;
			floatColor[0] = floatColor[1] = floatColor[2] = luminance;
		}
		else if (r_greyscale->value)
		{
			float luminance;

			luminance = LUMA(dl->color[0], dl->color[1], dl->color[2]) * 255.0f;
			floatColor[0] = LERP(dl->color[0] * 255.0f, luminance, r_greyscale->value);
			floatColor[1] = LERP(dl->color[1] * 255.0f, luminance, r_greyscale->value);
			floatColor[2] = LERP(dl->color[2] * 255.0f, luminance, r_greyscale->value);
		}
		else
		{
			floatColor[0] = dl->color[0] * 255.0f;
			floatColor[1] = dl->color[1] * 255.0f;
			floatColor[2] = dl->color[2] * 255.0f;
		}

		floatColorVec = _mm_loadu_ps(floatColor);
		//floatColorVec0 = vec_ld(0, floatColor);
		//floatColorVec1 = vec_ld(11, floatColor);
		//floatColorVec0 = vec_perm(floatColorVec0, floatColorVec0, floatColorVecPerm);
		for (i = 0; i < tess.numVertexes; i++, texCoords += 2/*, colors += 4*/) {
			__m128 colorVec0;
			__m128i colorInt;
			int		clip = 0;
			vec3_t	dist;

			VectorSubtract(origin, tess.xyz[i], dist);

			backEnd.pc.c_dlightVertexes++;

			texCoords[0] = 0.5f + dist[0] * scale;
			texCoords[1] = 0.5f + dist[1] * scale;

			if (!r_dlightBacks->integer &&
				// dist . tess.normal[i]
				(dist[0] * tess.normal[i][0] +
					dist[1] * tess.normal[i][1] +
					dist[2] * tess.normal[i][2]) < 0.0f) {
				clip = 63;
			}
			else {
				if (texCoords[0] < 0.0f) {
					clip |= 1;
				}
				else if (texCoords[0] > 1.0f) {
					clip |= 2;
				}
				if (texCoords[1] < 0.0f) {
					clip |= 4;
				}
				else if (texCoords[1] > 1.0f) {
					clip |= 8;
				}

				// modulate the strength based on the height and color
				if (dist[2] > radius) {
					clip |= 16;
					modulate = 0.0f;
				}
				else if (dist[2] < -radius) {
					clip |= 32;
					modulate = 0.0f;
				}
				else {
					dist[2] = Q_fabs(dist[2]);
					if (dist[2] < radius * 0.5f) {
						modulate = 1.0f;
					}
					else {
						modulate = 2.0f * (radius - dist[2]) * scale;
					}
				}
			}
			clipBits[i] = clip;

			modulateVec = _mm_set1_ps(modulate);
			modulateVec = _mm_add_ps(modulateVec, vSel1);
			//modulateVec = _mm_setr_ps(modulate, modulate, modulate, 255);
			//modulateVec = vec_ld(0, (float *)&modulate);
			//modulateVec = vec_perm(modulateVec, modulateVec, modulatePerm);
			colorVec0 = _mm_mul_ps(floatColorVec, modulateVec);
			colorVec0 = _mm_min_ps(colorVec0, vSel);
			//colorVec = vec_madd(floatColorVec0, modulateVec, zero);
			colorInt = _mm_cvtps_epi32(colorVec0);
			colorVec = _mm_slli_si128(colorVec, 1);
			colorVec = _mm_add_epi8(colorVec, colorInt);
			step++;
			if (step == 4)
			{
				step = 0;
				colorVec = _mm_shuffle_epi8(colorVec, shuf);
				_mm_store_si128((__m128i *)colors, colorVec);
				colorVec = _mm_setzero_si128();
				colors += 16;
			}
			//colorInt = vec_cts(colorVec, 0);	// RGBx
			//colorShort = vec_pack(colorInt, colorInt);		// RGBxRGBx
			//colorChar = vec_packsu(colorShort, colorShort);	// RGBxRGBxRGBxRGBx
			//colorChar = vec_sel(colorChar, vSel, vSel);		// RGBARGBARGBARGBA replace alpha with 255
			//vec_ste((vector unsigned int)colorChar, 0, (unsigned int *)colors);	// store color
		}
		if (step != 0)
		{
			colorVec = _mm_shuffle_epi8(colorVec, shuf);
			_mm_store_si128((__m128i *)colors, colorVec);
		}

		// build a list of triangles that need light
		numIndexes = 0;
		for (i = 0; i < tess.numIndexes; i += 3) {
			glIndex_t		a, b, c;

			a = tess.indexes[i];
			b = tess.indexes[i + 1];
			c = tess.indexes[i + 2];
			if (clipBits[a] & clipBits[b] & clipBits[c]) {
				continue;	// not lighted
			}
			hitIndexes[numIndexes] = a;
			hitIndexes[numIndexes + 1] = b;
			hitIndexes[numIndexes + 2] = c;
			numIndexes += 3;
		}

		if (!numIndexes) {
			continue;
		}
#if 0
		qglEnableClientState(GL_TEXTURE_COORD_ARRAY);
		qglTexCoordPointer(2, GL_FLOAT, 0, texCoordsArray[0]);

		qglEnableClientState(GL_COLOR_ARRAY);
		qglColorPointer(4, GL_UNSIGNED_BYTE, 0, colorArray);

		GL_Bind(tr.dlightImage);
		// include GLS_DEPTHFUNC_EQUAL so alpha tested surfaces don't add light
		// where they aren't rendered
		if (dl->additive) {
			GL_State(GLS_SRCBLEND_ONE | GLS_DSTBLEND_ONE | GLS_DEPTHFUNC_EQUAL);
		}
		else {
			GL_State(GLS_SRCBLEND_DST_COLOR | GLS_DSTBLEND_ONE | GLS_DEPTHFUNC_EQUAL);
		}
		R_DrawElements(numIndexes, hitIndexes);
#endif
		backEnd.pc.c_totalIndexes += numIndexes;
		backEnd.pc.c_dlightIndexes += numIndexes;
	}
}

struct dlight_s mydligths[100000];

#define TEST_SAMPLESIZE 140000
void maintest_projectdlighttexture(void)
{
	int i, j, k;
	Timer timer0, timer1, timer2;
	Cputime cpu0, cpu1, cpu2;
	double elapsed0, elapsed1, elapsed2;
	double cputime0, cputime1, cputime2;

	tess.dlightBits = 3;
	tess.numVertexes = 1;
	backEnd.refdef.num_dlights = 1;
	backEnd.refdef.dlights = mydligths;


	bdmpx_create(NULL, "bdmpx_dlighttexture.bin", BDMPX_OP_READ);

	for (j = 0; j < TEST_SAMPLESIZE; j++)
	{
		int ndl = 5;
		int nv = 5;
		int ni = 5;
		int dls = sizeof(mydligths);
		int xyz = sizeof(tess.xyz);
		int ids = sizeof(tess.indexes);
		int ercd = bdmpx_read(NULL, 6, &ndl, &backEnd.refdef.num_dlights, &dls, backEnd.refdef.dlights,
			&nv, &tess.numVertexes, &xyz, tess.xyz,
			&ni, &tess.numIndexes, &ids, tess.indexes);

		if (ercd == 0)
			break;
		assert(ercd == 6);
		assert(ndl == 4);
		assert(dls == backEnd.refdef.num_dlights*sizeof(struct dlight_s));
		assert(nv == 4);
		assert(ni == 4);
		assert(xyz == tess.numVertexes*sizeof(vec4_t));
		assert(ids == tess.numIndexes*sizeof(glIndex_t));


		cpu0.reset();
		timer0.reset();
		for (i = 0; i < 10; i++)
		{
			ProjectDlightTexture_scalar();
		}
		elapsed0 = timer0.accum();
		cputime0 = cpu0.accum();

		cpu1.reset();
		timer1.reset();
		for (i = 0; i < 10; i++)
		{
			ProjectDlightTexture_vector();
		}
		elapsed1 = timer1.accum();
		cputime1 = cpu1.accum();

	}

	printf("size %d\n", j);
	printf("scalar %4.4f %4.4f %4.4f %4.4f\n", elapsed0, cputime0, (elapsed1 - elapsed0) / elapsed1, (cputime1 - cputime0) / cputime1);
	printf("vector %4.4f %4.4f %4.4f %4.4f\n", elapsed1, cputime1, (elapsed0 - elapsed1) / elapsed0, (cputime0 - cputime1) / cputime0);
}