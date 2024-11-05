
#include <math.h>
#include <float.h>
#include <cassert>
#include <intrin.h>
#include <stdio.h>
#include <stdint.h>

#include "bdmpx.h"
#include "timing.h"
#include "platform.h"
#include "csv.h"

#define	myftol(x) ((int)(x))

#define DotProduct(x,y)			((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2])
#define VectorSubtract(a,b,c)	((c)[0]=(a)[0]-(b)[0],(c)[1]=(a)[1]-(b)[1],(c)[2]=(a)[2]-(b)[2])
#define VectorAdd(a,b,c)		((c)[0]=(a)[0]+(b)[0],(c)[1]=(a)[1]+(b)[1],(c)[2]=(a)[2]+(b)[2])
#define VectorCopy(a,b)			((b)[0]=(a)[0],(b)[1]=(a)[1],(b)[2]=(a)[2])
#define	VectorScale(v, s, o)	((o)[0]=(v)[0]*(s),(o)[1]=(v)[1]*(s),(o)[2]=(v)[2]*(s))
#define	VectorMA(v, s, b, o)	((o)[0]=(v)[0]+(b)[0]*(s),(o)[1]=(v)[1]+(b)[1]*(s),(o)[2]=(v)[2]+(b)[2]*(s))

typedef int		qhandle_t;

typedef float vec_t;
typedef vec_t vec2_t[2];
typedef vec_t vec3_t[3];
typedef vec_t vec4_t[4];

typedef enum { qfalse = 0, qtrue } qboolean;

typedef union floatint_u
{
	int i;
	unsigned int u;
	float f;
	byte b[4];
}
floatint_t;

typedef enum {
	RT_MODEL,
	RT_POLY,
	RT_SPRITE,
	RT_BEAM,
	RT_RAIL_CORE,
	RT_RAIL_RINGS,
	RT_LIGHTNING,
	RT_PORTALSURFACE,		// doesn't draw anything, just info for portals

	RT_MAX_REF_ENTITY_TYPE
} refEntityType_t;

typedef struct {
	refEntityType_t	reType;
	int			renderfx;

	qhandle_t	hModel;				// opaque type outside refresh

	// most recent data
	vec3_t		lightingOrigin;		// so multi-part models can be lit identically (RF_LIGHTING_ORIGIN)
	float		shadowPlane;		// projection shadows go here, stencils go slightly lower

	vec3_t		axis[3];			// rotation vectors
	qboolean	nonNormalizedAxes;	// axis are not normalized, i.e. they have scale
	float		origin[3];			// also used as MODEL_BEAM's "from"
	int			frame;				// also used as MODEL_BEAM's diameter

	// previous data for frame interpolation
	float		oldorigin[3];		// also used as MODEL_BEAM's "to"
	int			oldframe;
	float		backlerp;			// 0.0 = current, 1.0 = old

	// texturing
	int			skinNum;			// inline skin index
	qhandle_t	customSkin;			// NULL for default skin
	qhandle_t	customShader;		// use one image for the entire thing

	// misc
	byte		shaderRGBA[4];		// colors used by rgbgen entity shaders
	float		shaderTexCoord[2];	// texture coordinates used by tcMod entity modifiers

	// subtracted from refdef time to control effect start times
	floatint_t	shaderTime;

	// extra sprite information
	float		radius;
	float		rotation;
} refEntity_t;

// a trRefEntity_t has all the information passed in by
// the client game, as well as some locally derived info
typedef struct {
	refEntity_t	e;

	float		axisLength;		// compensate for non-normalized axis

	qboolean	needDlights;	// true for bmodels that touch a dlight
	qboolean	lightingCalculated;
	vec3_t		lightDir;		// normalized direction towards light
	vec3_t		ambientLight;	// color normalized to 0-255
	int			ambientLightInt;	// 32 bit rgba packed
	vec3_t		directedLight;
} trRefEntity_t;

#define	SHADER_MAX_VERTEXES	1000

trRefEntity_t entity;
static struct {
	trRefEntity_t	*currentEntity;
} backEnd = { &entity };

static struct {
	int			numVertexes;
	QALIGNA(16) vec4_t		normal[SHADER_MAX_VERTEXES];
} tess;

static void RB_CalcDiffuseColor( unsigned char *colors )
{
	int				i, j;
	float			/**v,*/ *normal;
	float			incoming;
	trRefEntity_t	*ent;
	int				ambientLightInt;
	vec3_t			ambientLight;
	vec3_t			lightDir;
	vec3_t			directedLight;
	int				numVertexes;
	ent = backEnd.currentEntity;
	ambientLightInt = ent->ambientLightInt;
	VectorCopy( ent->ambientLight, ambientLight );
	VectorCopy( ent->directedLight, directedLight );
	VectorCopy( ent->lightDir, lightDir );

	//v = tess.xyz[0];
	normal = tess.normal[0];

	numVertexes = tess.numVertexes;
	for (i = 0 ; i < numVertexes ; i++, /*v += 4,*/ normal += 4) {
		incoming = DotProduct (normal, lightDir);
		if ( incoming <= 0 ) {
			*(int *)&colors[i*4] = ambientLightInt;
			continue;
		} 
		j = myftol( ambientLight[0] + incoming * directedLight[0] );
		if ( j > 255 ) {
			j = 255;
		}
		colors[i*4+0] = j;

		j = myftol( ambientLight[1] + incoming * directedLight[1] );
		if ( j > 255 ) {
			j = 255;
		}
		colors[i*4+1] = j;

		j = myftol( ambientLight[2] + incoming * directedLight[2] );
		if ( j > 255 ) {
			j = 255;
		}
		colors[i*4+2] = j;

		colors[i*4+3] = 255;
	}
}

static void RB_CalcDiffuseColor_enhanced(unsigned char *colors)
{
	int				i, j;
	float			/**v,*/ *normal;
	float			incoming;
	trRefEntity_t	*ent;
	int				ambientLightInt;
	vec3_t			ambientLight;
	vec3_t			directedLight;
	int				numVertexes;
	ent = backEnd.currentEntity;
	ambientLightInt = ent->ambientLightInt;
	VectorCopy(ent->ambientLight, ambientLight);
	VectorCopy(ent->directedLight, directedLight);
	//VectorCopy(ent->lightDir, lightDir);
	__m128 normalVec0, lightDirVec, zero, incomingVec0;
	lightDirVec = _mm_loadu_ps(ent->lightDir);
	zero = _mm_setzero_ps();

	//v = tess.xyz[0];
	normal = tess.normal[0];

	numVertexes = tess.numVertexes;
	for (i = 0; i < numVertexes; i++, /*v += 4,*/ normal += 4) {
		normalVec0 = _mm_load_ps(normal);
#if 1
		incomingVec0 = _mm_dp_ps(normalVec0, lightDirVec, 0x77);
#else
		incomingVec0 = _mm_mul_ps(normalVec0, lightDirVec);
		incomingVec0 = _mm_hadd_ps(incomingVec0, zero);
		incomingVec0 = _mm_hadd_ps(incomingVec0, zero);
#endif
		//incoming = DotProduct(normal, lightDir);
		//if (incoming <= 0) {
		if(_mm_comile_ss(incomingVec0,zero)) {
			*(int *)&colors[i * 4] = ambientLightInt;
			continue;
		}
	    _mm_store_ss(&incoming, incomingVec0);
		j = myftol(ambientLight[0] + incoming * directedLight[0]);
		if (j > 255) {
			j = 255;
		}
		colors[i * 4 + 0] = j;

		j = myftol(ambientLight[1] + incoming * directedLight[1]);
		if (j > 255) {
			j = 255;
		}
		colors[i * 4 + 1] = j;

		j = myftol(ambientLight[2] + incoming * directedLight[2]);
		if (j > 255) {
			j = 255;
		}
		colors[i * 4 + 2] = j;

		colors[i * 4 + 3] = 255;
	}
}

/*
** RB_CalcDiffuseColor
**
** The basic vertex lighting calc
*/
void RB_CalcDiffuseColor_vector( unsigned char *colors )
{
  int				i;
  float			/**v,*/ *normal;
  trRefEntity_t	*ent;
  int				numVertexes, cycles;
  __m128 vSel2;
  __m128 ambientLightVec;
  __m128 directedLightVec;
  __m128 lightDirVec;
  __m128 normalVec0, normalVec1, normalVec2, normalVec3;
  __m128 incomingVec0, incomingVec1, incomingVec2, incomingVec3;
  __m128 zero, jVec0, jVec1, jVec2, jVec3;
  __m128i jVecInt0, jVecInt1, jVecInt2, jVecInt3;
  __m128i shuf0;
  /*
  vector unsigned char vSel = (vector unsigned char)(0x00, 0x00, 0x00, 0xff,
  0x00, 0x00, 0x00, 0xff,
  0x00, 0x00, 0x00, 0xff,
  0x00, 0x00, 0x00, 0xff);
  vector float ambientLightVec;
  vector float directedLightVec;
  vector float lightDirVec;
  vector float normalVec0, normalVec1;
  vector float incomingVec0, incomingVec1, incomingVec2;
  vector float zero, jVec;
  vector signed int jVecInt;
  vector signed short jVecShort;
  vector unsigned char jVecChar, normalPerm;*/
  ent = backEnd.currentEntity;
  // A lot of this could be simplified if we made sure
  // entities light info was 16-byte aligned.
  ambientLightVec = _mm_loadu_ps(ent->ambientLight);
  ambientLightVec = _mm_mul_ps(ambientLightVec, _mm_setr_ps(1.,1.,1.,0.));
  ambientLightVec = _mm_add_ps(ambientLightVec, _mm_setr_ps(0., 0., 0., 255.));
  /*jVecChar = vec_lvsl(0, ent->ambientLight);
  ambientLightVec = vec_ld(0, (vector float *)ent->ambientLight);
  jVec = vec_ld(11, (vector float *)ent->ambientLight);
  ambientLightVec = vec_perm(ambientLightVec,jVec,jVecChar);*/

  directedLightVec = _mm_loadu_ps(ent->directedLight);
  /*jVecChar = vec_lvsl(0, ent->directedLight);
  directedLightVec = vec_ld(0,(vector float *)ent->directedLight);
  jVec = vec_ld(11,(vector float *)ent->directedLight);
  directedLightVec = vec_perm(directedLightVec,jVec,jVecChar);*/

  lightDirVec = _mm_loadu_ps(ent->lightDir);
  /*jVecChar = vec_lvsl(0, ent->lightDir);
  lightDirVec = vec_ld(0,(vector float *)ent->lightDir);
  jVec = vec_ld(11,(vector float *)ent->lightDir);
  lightDirVec = vec_perm(lightDirVec,jVec,jVecChar);*/ 

  zero = _mm_setzero_ps();
  //zero4 = _mm_castsi128_ps(_mm_set_epi32(0, 0xffffffff, 0xffffffff, 0xffffffff));
  //maxu8 = _mm_set1_epi32(255);
  //vSel = _mm_set_ps(255, 0, 0, 0);
  vSel2 = _mm_set_ps(255, 255, 255, 255);
  //vCompact = _mm_set_epi8(0xff, 0xff, 0xff, 0xff,
  //  0xff, 0xff, 0xff, 0xff,
  //  0xff, 0xff, 0xff, 0xff,
  //  12, 8, 4, 0);
  //vColorMask = _mm_set_epi8(0, 0, 0, 0,
  //  0, 0, 0, 0,
  //  0, 0, 0, 0,
  //  0xff, 0xff, 0xff, 0xff);
  //zero = (vector float)vec_splat_s8(0);

  //v = tess.xyz[0];
  normal = tess.normal[0];

  //normalPerm = vec_lvsl(0,normal);
  numVertexes = tess.numVertexes;
  //normalVec = _mm_loadu_ps(normal);

  //bdmpx_write(NULL, 4,sizeof(vec3_t),ent->ambientLight,sizeof(vec3_t),ent->directedLight,sizeof(vec3_t),ent->lightDir,sizeof(vec4_t)*numVertexes,normal);
  //assert((normal)[3] == 0);
  //normalVec = _mm_and_ps(zero4, normalVec);
  for (i = 0 ; i < numVertexes; i+=4, /*v += 4,*/ normal += 16)
  {
    normalVec0 = _mm_load_ps(normal);
    /*normalVec0 = vec_ld(0,(vector float *)normal);
    normalVec1 = vec_ld(11,(vector float *)normal);
    normalVec0 = vec_perm(normalVec0,normalVec1,normalPerm);*/
    incomingVec0 = _mm_mul_ps(normalVec0, lightDirVec);
    //incomingVec0 = vec_madd(normalVec0, lightDirVec, zero);
    incomingVec0 = _mm_hadd_ps(incomingVec0, zero);
    incomingVec0 = _mm_hadd_ps(incomingVec0, zero);
    //incomingVec1 = vec_sld(incomingVec0,incomingVec0,4);
    //incomingVec2 = vec_add(incomingVec0,incomingVec1);
    //incomingVec1 = vec_sld(incomingVec1,incomingVec1,4);
    //incomingVec2 = vec_add(incomingVec2,incomingVec1);
    incomingVec0 = _mm_shuffle_ps(incomingVec0, incomingVec0, 0xc0);
    //incomingVec0 = vec_splat(incomingVec2,0);
    incomingVec0 = _mm_max_ps(incomingVec0,zero);
    //incomingVec0 = vec_max(incomingVec0,zero);
    //assert((normal+4)[3] == 0);
    //normalVec = _mm_and_ps(zero4, normalVec);
    //normalPerm = vec_lvsl(12,normal);
    jVec0 = _mm_mul_ps(incomingVec0, directedLightVec);
    jVec0 = _mm_add_ps(jVec0, ambientLightVec);
    //jVec = vec_madd(incomingVec0, directedLightVec, ambientLightVec);
    jVec0 = _mm_min_ps(jVec0, vSel2);

    jVecInt0 = _mm_cvttps_epi32(jVec0);

    normalVec1 = _mm_load_ps(normal+4);
    incomingVec1 = _mm_mul_ps(normalVec1, lightDirVec);
    incomingVec1 = _mm_hadd_ps(incomingVec1, zero);
    incomingVec1 = _mm_hadd_ps(incomingVec1, zero);
    incomingVec1 = _mm_shuffle_ps(incomingVec1, incomingVec1, 0xc0);
    incomingVec1 = _mm_max_ps(incomingVec1,zero);
    jVec1 = _mm_mul_ps(incomingVec1, directedLightVec);
    jVec1 = _mm_add_ps(jVec1, ambientLightVec);
    jVec1 = _mm_min_ps(jVec1, vSel2);
    jVecInt1 = _mm_cvttps_epi32(jVec1);
    jVecInt1 = _mm_slli_si128(jVecInt1,1);

    jVecInt0 = _mm_add_epi8(jVecInt0, jVecInt1);
    
    normalVec2 = _mm_load_ps(normal+8);
    incomingVec2 = _mm_mul_ps(normalVec2, lightDirVec);
    incomingVec2 = _mm_hadd_ps(incomingVec2, zero);
    incomingVec2 = _mm_hadd_ps(incomingVec2, zero);
    incomingVec2 = _mm_shuffle_ps(incomingVec2, incomingVec2, 0xc0);
    incomingVec2 = _mm_max_ps(incomingVec2,zero);
    jVec2 = _mm_mul_ps(incomingVec2, directedLightVec);
    jVec2 = _mm_add_ps(jVec2, ambientLightVec);
    jVec2 = _mm_min_ps(jVec2, vSel2);
    jVecInt2 = _mm_cvttps_epi32(jVec2);
    jVecInt2 = _mm_slli_si128(jVecInt2,2);

    jVecInt0 = _mm_add_epi8(jVecInt0, jVecInt2);
    
    normalVec3 = _mm_load_ps(normal+12);
    incomingVec3 = _mm_mul_ps(normalVec3, lightDirVec);
    incomingVec3 = _mm_hadd_ps(incomingVec3, zero);
    incomingVec3 = _mm_hadd_ps(incomingVec3, zero);
    incomingVec3 = _mm_shuffle_ps(incomingVec3, incomingVec3, 0xc0);
    incomingVec3 = _mm_max_ps(incomingVec3,zero);
    jVec3 = _mm_mul_ps(incomingVec3, directedLightVec);
    jVec3 = _mm_add_ps(jVec3, ambientLightVec);
    jVec3 = _mm_min_ps(jVec3, vSel2);
    jVecInt3 = _mm_cvttps_epi32(jVec3);
    jVecInt3 = _mm_slli_si128(jVecInt3,3);
    
    jVecInt0 = _mm_add_epi8(jVecInt0,jVecInt3);
    shuf0 = _mm_setr_epi8(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15);
    jVecInt0 = _mm_shuffle_epi8(jVecInt0,shuf0);
    //colors[i * 4] = _mm_cvt_ss2si(jVec0);
    //colors[i * 4 + 1] = _mm_cvt_ss2si(_mm_movehdup_ps(jVec0));
    //colors[i * 4 + 2] = _mm_cvt_ss2si(_mm_unpackhi_ps(jVec0,jVec0));
    //colors[i * 4 + 3] = 255;
    //jVecInt = _mm_cvttps_epi32(jVec);
    //jVecInt = vec_cts(jVec,0);	// RGBx
    //jVecCmp = _mm_cmplt_epi32(jVecInt, vSel);
    //jVecInt = _mm_and_si128(jVecCmp, jVecInt);
    //jVecInt = _mm_add_epi32(jVecInt, _mm_andnot_si128(jVecCmp, vSel2));
    //jVecInt = _mm_shuffle_epi8(jVecInt, vCompact);
    //jVecShort = vec_pack(jVecInt,jVecInt);		// RGBxRGBx
    //jVecChar = vec_packsu(jVecShort,jVecShort);	// RGBxRGBxRGBxRGBx
    //jVecChar = vec_sel(jVecChar,vSel,vSel);		// RGBARGBARGBARGBA replace alpha with 255
    //_mm_maskmoveu_si128(jVecInt, vColorMask, (char*)&colors[i*4]);
    _mm_storeu_si128((__m128i *)&colors[i * 4], jVecInt0);
    //vec_ste((vector unsigned int)jVecChar,0,(unsigned int *)&colors[i*4]);	// store color
  }
  //_mm_sfence();
}

void RB_CalcDiffuseColor_sse4(unsigned char *colors)
{
	int				i;
	float			/**v,*/ *normal;
	trRefEntity_t	*ent;
	int				numVertexes, cycles;
	__m128 vSel2;
	__m128 ambientLightVec;
	__m128 directedLightVec;
	__m128 lightDirVec;
	__m128 normalVec0, normalVec1, normalVec2, normalVec3;
	__m128 incomingVec0, incomingVec1, incomingVec2, incomingVec3;
	__m128 zero, jVec0, jVec1, jVec2, jVec3;
	__m128i jVecInt0, jVecInt1, jVecInt2, jVecInt3;
	__m128i shuf0;
	
	ent = backEnd.currentEntity;
	ambientLightVec = _mm_loadu_ps(ent->ambientLight);
	ambientLightVec = _mm_mul_ps(ambientLightVec, _mm_setr_ps(1., 1., 1., 0.));
	ambientLightVec = _mm_add_ps(ambientLightVec, _mm_setr_ps(0., 0., 0., 255.));

	directedLightVec = _mm_loadu_ps(ent->directedLight);

	lightDirVec = _mm_loadu_ps(ent->lightDir);

	zero = _mm_setzero_ps();
	vSel2 = _mm_set_ps(255, 255, 255, 255);

	//v = tess.xyz[0];
	normal = tess.normal[0];

	numVertexes = tess.numVertexes;
	//normalVec = _mm_loadu_ps(normal);

	//bdmpx_write(NULL, 4,sizeof(vec3_t),ent->ambientLight,sizeof(vec3_t),ent->directedLight,sizeof(vec3_t),ent->lightDir,sizeof(vec4_t)*numVertexes,normal);
	//assert((normal)[3] == 0);
	//normalVec = _mm_and_ps(zero4, normalVec);
	for (i = 0; i < numVertexes; i += 4, /*v += 4,*/ normal += 16)
	{
		normalVec0 = _mm_load_ps(normal);
		incomingVec0 = _mm_dp_ps(normalVec0, lightDirVec, 0x77);
		incomingVec0 = _mm_max_ps(incomingVec0, zero);
		jVec0 = _mm_mul_ps(incomingVec0, directedLightVec);
		jVec0 = _mm_add_ps(jVec0, ambientLightVec);
		jVec0 = _mm_min_ps(jVec0, vSel2);

		jVecInt0 = _mm_cvttps_epi32(jVec0);

		normalVec1 = _mm_load_ps(normal + 4);
		incomingVec1 = _mm_dp_ps(normalVec1, lightDirVec, 0x77);
		incomingVec1 = _mm_max_ps(incomingVec1, zero);
		jVec1 = _mm_mul_ps(incomingVec1, directedLightVec);
		jVec1 = _mm_add_ps(jVec1, ambientLightVec);
		jVec1 = _mm_min_ps(jVec1, vSel2);
		jVecInt1 = _mm_cvttps_epi32(jVec1);
		jVecInt1 = _mm_slli_si128(jVecInt1, 1);

		jVecInt0 = _mm_add_epi8(jVecInt0, jVecInt1);

		normalVec2 = _mm_load_ps(normal + 8);
		incomingVec2 = _mm_dp_ps(normalVec2, lightDirVec, 0x77);
		incomingVec2 = _mm_max_ps(incomingVec2, zero);
		jVec2 = _mm_mul_ps(incomingVec2, directedLightVec);
		jVec2 = _mm_add_ps(jVec2, ambientLightVec);
		jVec2 = _mm_min_ps(jVec2, vSel2);
		jVecInt2 = _mm_cvttps_epi32(jVec2);
		jVecInt2 = _mm_slli_si128(jVecInt2, 2);

		jVecInt0 = _mm_add_epi8(jVecInt0, jVecInt2);

		normalVec3 = _mm_load_ps(normal + 12);
		incomingVec3 = _mm_dp_ps(normalVec3, lightDirVec, 0x77);
		incomingVec3 = _mm_max_ps(incomingVec3, zero);
		jVec3 = _mm_mul_ps(incomingVec3, directedLightVec);
		jVec3 = _mm_add_ps(jVec3, ambientLightVec);
		jVec3 = _mm_min_ps(jVec3, vSel2);
		jVecInt3 = _mm_cvttps_epi32(jVec3);
		jVecInt3 = _mm_slli_si128(jVecInt3, 3);

		jVecInt0 = _mm_add_epi8(jVecInt0, jVecInt3);
		shuf0 = _mm_setr_epi8(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15);
		jVecInt0 = _mm_shuffle_epi8(jVecInt0, shuf0);
		_mm_storeu_si128((__m128i *)&colors[i * 4], jVecInt0);
	}
}

typedef void (*diffuse_fn)(unsigned char* outputs);

struct function_data
{
	void* fp;
	const char* fname;
};

static struct function_data data_fn[] =
{
	{ RB_CalcDiffuseColor,            "    diffuse" },
#if idsse
	{ RB_CalcDiffuseColor_enhanced,   "diffuse_enh" },
	{ RB_CalcDiffuseColor_vector,     "diffuse_vec" },
	{ RB_CalcDiffuseColor_sse4,       " diffuse_dp" },
#endif
	{ 0, 0 }
};

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(X[0]))

#define MAX_COLOR_SAMPLES SHADER_MAX_VERTEXES
struct {
	QALIGNA(16)  uint32_t thecolors[(ARRAY_SIZE(data_fn) - 1) * MAX_COLOR_SAMPLES];
} data;

#define TEST_SAMPLESIZE 100000

#define MAX_ERRORS 20

#define REPETITIONS 1
void maintest_diffusecolor(void)
{
  int sz, ercd;
  int i, j, k;
  int err = 0;
  int tested = 0;
  Timer timer[ARRAY_SIZE(data_fn) - 1];
  double elapsed[ARRAY_SIZE(data_fn) - 1];

  ercd = csv_open("./results_diffusecolor.csv");
  if (ercd != 0)
  {
	  printf("Could not open the csv file.\n\n");
  }

  //bdmpx_set_option(BDMPX_OPTION_PRECACHE_FILE_FOR_READ, 1);

  ercd = bdmpx_create(NULL, "bdmpx_diffusecolor.bin", BDMPX_OP_READ);
  if (ercd != 0)
  {
	  printf("Could not open the test input file. Aborting.\n");
	  return;
  }

  for(j = 0; j < TEST_SAMPLESIZE ; j++)
  {
    //accum.reset();
    sz = sizeof(tess.normal);
    ercd = bdmpx_read(NULL, 4, NULL, entity.ambientLight, NULL, entity.directedLight, NULL, entity.lightDir, &sz, tess.normal);
	tess.numVertexes = sz/sizeof(vec4_t);

	entity.ambientLightInt = ((int)entity.ambientLight[0]) +
		((int)entity.ambientLight[1] << 8) + ((int)entity.ambientLight[2] << 16) + (255 << 24);

    if(tess.numVertexes > MAX_COLOR_SAMPLES)
    {
      printf("Too many tess vertexes: %d %d", j , tess.numVertexes);
	  return;
    }

    if(ercd != 4)
      break;

	uint32_t* out;

	memset(&data.thecolors, 0, sizeof(data.thecolors));

	for (tested = 0, out = data.thecolors; data_fn[tested].fp != 0; tested++, out += MAX_COLOR_SAMPLES)
	{
		diffuse_fn callme = (diffuse_fn)data_fn[tested].fp;
		const char* myinfo = data_fn[tested].fname;

		timer[tested].reset();
		for (int r = 0; r < REPETITIONS; r++)
			callme((unsigned char*)out);
		elapsed[tested] = timer[tested].accum();
	}

	for (i = 1; i < tested; i++)
	{
		for (k = 0; k < tess.numVertexes/*MAX_COLOR_SAMPLES*/ && err < MAX_ERRORS; k++)
		{
			if (data.thecolors[k] != data.thecolors[i* MAX_COLOR_SAMPLES + k])
			{
				printf("id:%d seg:%d sam:%d/%d orig:%x here:%x\n", i, j, k, tess.numVertexes, data.thecolors[k], data.thecolors[i * MAX_COLOR_SAMPLES + k]);
				err++;
			}
		}
	}
  }
  
  printf("Diffusecolor Test segments %d\n", j);
  for (k = 0; k < tested; k++)
  {
	  const char* myinfo = data_fn[k].fname;
	  printf("%d %s: %4.4f %4.4f\n", k, myinfo, elapsed[k], elapsed[k] / elapsed[0]);
	  csv_put_float(elapsed[k]);
  }

  csv_put_string(",");
  for (k = 1; k < tested; k++)
  {
	  csv_put_float(elapsed[k] / elapsed[0]);
  }

  csv_put_string("\n");
  csv_close();

  printf("\n");

}