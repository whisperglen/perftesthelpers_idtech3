
#include <math.h>
#include <float.h>
#include <cassert>
#include <intrin.h>
#include <stdio.h>

#include "bdmpx.h"
#include "timing.h"
#include "platform.h"
#include "csv.h"

#define	MAX_QPATH			64		// max length of a quake game pathname
#define	SHADER_MAX_VERTEXES	1000
// vertex scales
#define	MD3_XYZ_SCALE		(1.0f/64)
#define	FOG_TABLE_SIZE		256
#define FUNCTABLE_SIZE		1024
#define FUNCTABLE_SIZE2		10
#define FUNCTABLE_MASK		(FUNCTABLE_SIZE-1)

#define ID_INLINE __inline

#ifndef M_PI
#define M_PI		3.14159265358979323846f	// matches value in gcc v2 math.h
#endif
#define DEG2RAD( a ) ( ( (a) * M_PI ) / 180.0F )

typedef float vec_t;
typedef vec_t vec2_t[2];
typedef vec_t vec3_t[3];
typedef vec_t vec4_t[4];
typedef vec_t vec5_t[5];

typedef	int	fixed4_t;
typedef	int	fixed8_t;
typedef	int	fixed16_t;

/*
** md3Surface_t
**
** CHUNK			SIZE
** header			sizeof( md3Surface_t )
** shaders			sizeof( md3Shader_t ) * numShaders
** triangles[0]		sizeof( md3Triangle_t ) * numTriangles
** st				sizeof( md3St_t ) * numVerts
** XyzNormals		sizeof( md3XyzNormal_t ) * numVerts * numFrames
*/
typedef struct {
	int		numVerts;
} md3Surface_t;

static struct shaderCommands_s
{
	int			numVertexes;
	QALIGNA(16) vec4_t		xyz[2 * SHADER_MAX_VERTEXES];
	QALIGNA(16) vec4_t		normal[2 * SHADER_MAX_VERTEXES];
	unsigned int guard;
} tess;

static struct
{
  float					sinTable[FUNCTABLE_SIZE];
} tr;

#define DotProduct(x,y)			((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2])

static ID_INLINE float Q_rsqrt( float number )
{
  __m128 number_f;
  float ret;

  if(number == 0.f) return 1.f / FLT_MIN;

  //number_f = _mm_load_ss(&number);
  number_f = _mm_set_ss(number);

  number_f = _mm_rsqrt_ss(number_f);

  _mm_store_ss(&ret, number_f);
  return ret;
}

// fast vector normalize routine that does not check to make sure
// that length != 0, nor does it return length, uses rsqrt approximation
static ID_INLINE void VectorNormalizeFast( vec3_t v )
{
	float ilength;

	ilength = Q_rsqrt( DotProduct( v, v ) );

	v[0] *= ilength;
	v[1] *= ilength;
	v[2] *= ilength;
}

static void VectorArrayNormalize(vec4_t *normals, unsigned int count)
{
//    assert(count);
	// given the input, it's safe to call VectorNormalizeFast
    while ( count-- ) {
        VectorNormalizeFast(normals[0]);
        normals++;
    }
}


#define SURFZS 100000
short mynewxyz[SURFZS];
short myoldxyz[SURFZS];
/*
** LerpMeshVertexes
*/
static void LerpMeshVertexes(md3Surface_t *surf, float backlerp)
{
	short	*oldXyz, *newXyz, *oldNormals, *newNormals;
	float	*outXyz, *outNormal;
	float	oldXyzScale, newXyzScale;
	float	oldNormalScale, newNormalScale;
	int		vertNum;
	unsigned lat, lng;
	int		numVerts;

	outXyz = tess.xyz[tess.numVertexes];
	outNormal = tess.normal[tess.numVertexes];

	newXyz = &mynewxyz[0];//(short *)((byte *)surf + surf->ofsXyzNormals)
		//+ (backEnd.currentEntity->e.frame * surf->numVerts * 4);
	newNormals = newXyz + 3;

	newXyzScale = MD3_XYZ_SCALE * (1.0f - backlerp);
	newNormalScale = 1.0f - backlerp;

	numVerts = surf->numVerts;

  //bdmpx_write(NULL, 4,sizeof(float),&backlerp,numVerts*sizeof(vec4_t),outXyz,numVerts*sizeof(vec4_t),outNormal,numVerts*8,newXyz);

	if ( backlerp == 0 ) {
		//
		// just copy the vertexes
		//
		for (vertNum=0 ; vertNum < numVerts ; vertNum++,
			newXyz += 4, newNormals += 4,
			outXyz += 4, outNormal += 4) 
		{

			outXyz[0] = newXyz[0] * newXyzScale;
			outXyz[1] = newXyz[1] * newXyzScale;
			outXyz[2] = newXyz[2] * newXyzScale;

			lat = ( newNormals[0] >> 8 ) & 0xff;
			lng = ( newNormals[0] & 0xff );
			lat *= (FUNCTABLE_SIZE/256);
			lng *= (FUNCTABLE_SIZE/256);

			// decode X as cos( lat ) * sin( long )
			// decode Y as sin( lat ) * sin( long )
			// decode Z as cos( long )

			outNormal[0] = tr.sinTable[(lat+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK] * tr.sinTable[lng];
			outNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
			outNormal[2] = tr.sinTable[(lng+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK];
		}
	} else {
		//
		// interpolate and copy the vertex and normal
		//
		oldXyz = &myoldxyz[0];//(short *)((byte *)surf + surf->ofsXyzNormals)
			//+ (backEnd.currentEntity->e.oldframe * surf->numVerts * 4);
		oldNormals = oldXyz + 3;

		oldXyzScale = MD3_XYZ_SCALE * backlerp;
		oldNormalScale = backlerp;

		//bdmpx_write(NULL, 1,numVerts*8,oldXyz);

		for (vertNum=0 ; vertNum < numVerts ; vertNum++,
			oldXyz += 4, newXyz += 4, oldNormals += 4, newNormals += 4,
			outXyz += 4, outNormal += 4) 
		{
			vec3_t uncompressedOldNormal, uncompressedNewNormal;

			// interpolate the xyz
			outXyz[0] = oldXyz[0] * oldXyzScale + newXyz[0] * newXyzScale;
			outXyz[1] = oldXyz[1] * oldXyzScale + newXyz[1] * newXyzScale;
			outXyz[2] = oldXyz[2] * oldXyzScale + newXyz[2] * newXyzScale;

			// FIXME: interpolate lat/long instead?
			lat = ( newNormals[0] >> 8 ) & 0xff;
			lng = ( newNormals[0] & 0xff );
			lat *= 4;
			lng *= 4;
			uncompressedNewNormal[0] = tr.sinTable[(lat+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK] * tr.sinTable[lng];
			uncompressedNewNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
			uncompressedNewNormal[2] = tr.sinTable[(lng+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK];

			lat = ( oldNormals[0] >> 8 ) & 0xff;
			lng = ( oldNormals[0] & 0xff );
			lat *= 4;
			lng *= 4;

			uncompressedOldNormal[0] = tr.sinTable[(lat+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK] * tr.sinTable[lng];
			uncompressedOldNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
			uncompressedOldNormal[2] = tr.sinTable[(lng+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK];

			outNormal[0] = uncompressedOldNormal[0] * oldNormalScale + uncompressedNewNormal[0] * newNormalScale;
			outNormal[1] = uncompressedOldNormal[1] * oldNormalScale + uncompressedNewNormal[1] * newNormalScale;
			outNormal[2] = uncompressedOldNormal[2] * oldNormalScale + uncompressedNewNormal[2] * newNormalScale;

//			VectorNormalize (outNormal);
		}
    	VectorArrayNormalize((vec4_t *)tess.normal[tess.numVertexes], numVerts);
   	}
}

/*
** LerpMeshVertexes
*/
#define idppc_altivec 1
static void LerpMeshVertexes_vector (md3Surface_t *surf, float backlerp) 
{
	short	*oldXyz, *newXyz, *oldNormals, *newNormals;
	float	*outXyz, *outNormal;
	float	oldXyzScale, newXyzScale;
	float	oldNormalScale, newNormalScale;
	int		vertNum;
	unsigned lat, lng;
	int		numVerts;

	outXyz = tess.xyz[tess.numVertexes];
	outNormal = tess.normal[tess.numVertexes];

	newXyz = &mynewxyz[0];//(short *)((byte *)surf + surf->ofsXyzNormals)
		//+ (backEnd.currentEntity->e.frame * surf->numVerts * 4);
	newNormals = newXyz + 3;

	newXyzScale = MD3_XYZ_SCALE * (1.0f - backlerp);
	newNormalScale = 1.0f - backlerp;

	numVerts = surf->numVerts;

	if ( backlerp == 0 ) {
#if idppc_altivec
    __m128 newXyzScaleVec;
    __m128i newNormalsIntVec, newNormalsIntVec0, newNormalsIntVec1;
    __m128 newNormalsFloatVec0, newNormalsFloatVec1;
    //const __m128i magicInt = _mm_set1_epi16(0x4B00);
    //const __m128 magicFloat = _mm_set1_ps(8388608.0f);
    //__m128 zero;

    newXyzScaleVec = _mm_setr_ps(newXyzScale, newXyzScale, newXyzScale, 0.f);
    //zero = _mm_setzero_ps();

		//vector signed short newNormalsVec0;
		//vector signed short newNormalsVec1;
		//vector signed int newNormalsIntVec;
		//vector float newNormalsFloatVec;
		//vector float newXyzScaleVec;
		//vector unsigned char newNormalsLoadPermute;
		//vector unsigned char newNormalsStorePermute;
		//vector float zero;
		
		//newNormalsStorePermute = vec_lvsl(0,(float *)&newXyzScaleVec);
		//newXyzScaleVec = *(vector float *)&newXyzScale;
		//newXyzScaleVec = vec_perm(newXyzScaleVec,newXyzScaleVec,newNormalsStorePermute);
		//newXyzScaleVec = vec_splat(newXyzScaleVec,0);		
		//newNormalsLoadPermute = vec_lvsl(0,newXyz);
		//newNormalsStorePermute = vec_lvsr(0,outXyz);
		//zero = (vector float)vec_splat_s8(0);
		//
		// just copy the vertexes
		//
		for (vertNum=0 ; vertNum < numVerts ; vertNum+=2,
			newXyz += 8, newNormals += 8,
			outXyz += 8, outNormal += 8) 
		{
      newNormalsIntVec = _mm_lddqu_si128((__m128i const*)newXyz);
      newNormalsIntVec0 = _mm_unpacklo_epi16(_mm_setzero_si128(), newNormalsIntVec);
      newNormalsIntVec0 = _mm_srai_epi32(newNormalsIntVec0, 16);
      newNormalsFloatVec0 = _mm_cvtepi32_ps(newNormalsIntVec0);
      //newNormalsIntVec0 = _mm_unpacklo_epi16(newNormalsIntVec, magicInt);
      //newNormalsFloatVec0 = _mm_sub_ps(_mm_castsi128_ps(newNormalsIntVec0), magicFloat);
			//newNormalsLoadPermute = vec_lvsl(0,newXyz);
			//newNormalsStorePermute = vec_lvsr(0,outXyz);
			//newNormalsVec0 = vec_ld(0,newXyz);
			//newNormalsVec1 = vec_ld(16,newXyz);
			//newNormalsVec0 = vec_perm(newNormalsVec0,newNormalsVec1,newNormalsLoadPermute);
			//newNormalsIntVec = vec_unpackh(newNormalsVec0);
			//newNormalsFloatVec = vec_ctf(newNormalsIntVec,0);
      newNormalsFloatVec0 = _mm_mul_ps(newNormalsFloatVec0, newXyzScaleVec);
			//newNormalsFloatVec = vec_madd(newNormalsFloatVec,newXyzScaleVec,zero);
			//newNormalsFloatVec = vec_perm(newNormalsFloatVec,newNormalsFloatVec,newNormalsStorePermute);
			//outXyz[0] = newXyz[0] * newXyzScale;
			//outXyz[1] = newXyz[1] * newXyzScale;
			//outXyz[2] = newXyz[2] * newXyzScale;

			lat = ( newNormals[0] >> 8 ) & 0xff;
			lng = ( newNormals[0] & 0xff );
			lat *= (FUNCTABLE_SIZE/256);
			lng *= (FUNCTABLE_SIZE/256);

			// decode X as cos( lat ) * sin( long )
			// decode Y as sin( lat ) * sin( long )
			// decode Z as cos( long )

			outNormal[0] = tr.sinTable[(lat+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK] * tr.sinTable[lng];
			outNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
			outNormal[2] = tr.sinTable[(lng+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK];

      newNormalsIntVec1 = _mm_unpackhi_epi16(_mm_setzero_si128(), newNormalsIntVec);
      newNormalsIntVec1 = _mm_srai_epi32(newNormalsIntVec1, 16);
      newNormalsFloatVec1 = _mm_cvtepi32_ps(newNormalsIntVec1);
      newNormalsFloatVec1 = _mm_mul_ps(newNormalsFloatVec1, newXyzScaleVec);

			lat = ( newNormals[4] >> 8 ) & 0xff;
			lng = ( newNormals[4] & 0xff );
			lat *= (FUNCTABLE_SIZE/256);
			lng *= (FUNCTABLE_SIZE/256);

			// decode X as cos( lat ) * sin( long )
			// decode Y as sin( lat ) * sin( long )
			// decode Z as cos( long )

			outNormal[0+4] = tr.sinTable[(lat+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK] * tr.sinTable[lng];
			outNormal[1+4] = tr.sinTable[lat] * tr.sinTable[lng];
			outNormal[2+4] = tr.sinTable[(lng+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK];

      _mm_storeu_ps(outXyz, newNormalsFloatVec0);
      _mm_storeu_ps(outXyz+4, newNormalsFloatVec1);
			//vec_ste(newNormalsFloatVec,0,outXyz);
			//vec_ste(newNormalsFloatVec,4,outXyz);
			//vec_ste(newNormalsFloatVec,8,outXyz);
		}
		
#else
		//
		// just copy the vertexes
		//
		for (vertNum=0 ; vertNum < numVerts ; vertNum++,
			newXyz += 4, newNormals += 4,
			outXyz += 4, outNormal += 4) 
		{

			outXyz[0] = newXyz[0] * newXyzScale;
			outXyz[1] = newXyz[1] * newXyzScale;
			outXyz[2] = newXyz[2] * newXyzScale;

			lat = ( newNormals[0] >> 8 ) & 0xff;
			lng = ( newNormals[0] & 0xff );
			lat *= (FUNCTABLE_SIZE/256);
			lng *= (FUNCTABLE_SIZE/256);

			// decode X as cos( lat ) * sin( long )
			// decode Y as sin( lat ) * sin( long )
			// decode Z as cos( long )

			outNormal[0] = tr.sinTable[(lat+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK] * tr.sinTable[lng];
			outNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
			outNormal[2] = tr.sinTable[(lng+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK];
		}
#endif
	} else {
		//
		// interpolate and copy the vertex and normal
		//
		oldXyz = &myoldxyz[0];//(short *)((byte *)surf + surf->ofsXyzNormals)
			//+ (backEnd.currentEntity->e.oldframe * surf->numVerts * 4);
		oldNormals = oldXyz + 3;

		oldXyzScale = MD3_XYZ_SCALE * backlerp;
		oldNormalScale = backlerp;

		for (vertNum=0 ; vertNum < numVerts ; vertNum++,
			oldXyz += 4, newXyz += 4, oldNormals += 4, newNormals += 4,
			outXyz += 4, outNormal += 4) 
		{
			vec3_t uncompressedOldNormal, uncompressedNewNormal;

			// interpolate the xyz
			outXyz[0] = oldXyz[0] * oldXyzScale + newXyz[0] * newXyzScale;
			outXyz[1] = oldXyz[1] * oldXyzScale + newXyz[1] * newXyzScale;
			outXyz[2] = oldXyz[2] * oldXyzScale + newXyz[2] * newXyzScale;

			// FIXME: interpolate lat/long instead?
			lat = ( newNormals[0] >> 8 ) & 0xff;
			lng = ( newNormals[0] & 0xff );
			lat *= 4;
			lng *= 4;
			uncompressedNewNormal[0] = tr.sinTable[(lat+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK] * tr.sinTable[lng];
			uncompressedNewNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
			uncompressedNewNormal[2] = tr.sinTable[(lng+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK];

			lat = ( oldNormals[0] >> 8 ) & 0xff;
			lng = ( oldNormals[0] & 0xff );
			lat *= 4;
			lng *= 4;

			uncompressedOldNormal[0] = tr.sinTable[(lat+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK] * tr.sinTable[lng];
			uncompressedOldNormal[1] = tr.sinTable[lat] * tr.sinTable[lng];
			uncompressedOldNormal[2] = tr.sinTable[(lng+(FUNCTABLE_SIZE/4))&FUNCTABLE_MASK];

			outNormal[0] = uncompressedOldNormal[0] * oldNormalScale + uncompressedNewNormal[0] * newNormalScale;
			outNormal[1] = uncompressedOldNormal[1] * oldNormalScale + uncompressedNewNormal[1] * newNormalScale;
			outNormal[2] = uncompressedOldNormal[2] * oldNormalScale + uncompressedNewNormal[2] * newNormalScale;

//			VectorNormalize (outNormal);
		}
    	VectorArrayNormalize((vec4_t *)tess.normal[tess.numVertexes], numVerts);
   	}
}

typedef void(*lerpmesh_fn)(void *surf, float backlerp);

struct function_data
{
	void* fp;
	const char* fname;
};

static struct function_data data_fn[] =
{
	{ LerpMeshVertexes,            "    lerpmesh" },
#if idsse
	{ LerpMeshVertexes_vector,     "lerpmesh_vec" },
#endif
	{ 0, 0 }
};

#define ARRAY_SIZE(X) (sizeof(X)/sizeof(X[0]))

#define TEST_SAMPLESIZE 140000

#define REPETITIONS 1
#define MAX_ERRORS 10
void maintest_lerpmeshvertexes(void)
{
  int i, j, k;
  int tested = 0;
  int ercd;
  int err = 0;
  Timer timer[ARRAY_SIZE(data_fn) - 1];
  Cputime cpu[ARRAY_SIZE(data_fn) - 1];
  double elapsed[ARRAY_SIZE(data_fn) - 1];
  double cputime[ARRAY_SIZE(data_fn) - 1];
  int backlerp_was_zero = 0;

  float backlerp;
  md3Surface_t surf;
  surf.numVerts = 0;

  for ( i = 0; i < FUNCTABLE_SIZE; i++ )
  {
		tr.sinTable[i]		= sin( DEG2RAD( i * 360.0f / ( ( float ) ( FUNCTABLE_SIZE - 1 ) ) ) );
  }

  if ((ARRAY_SIZE(data_fn) - 1) * SHADER_MAX_VERTEXES > min(ARRAY_SIZE(tess.xyz), ARRAY_SIZE(tess.normal)))
  {
	  printf("Increase tess.xyz and tess.normal to: %d * SHADER_MAX_VERTEXES\n", ARRAY_SIZE(data_fn) - 1);
	  return;
  }

  ercd = csv_open("./results_lerpmesh.csv");
  if (ercd != 0)
  {
	  printf("Could not open the csv file.\n\n");
  }

  ercd = bdmpx_create(NULL, "bdmpx_lerp.bin", BDMPX_OP_READ);
  if (ercd != 0)
  {
	  printf("Could not open the test input file. Aborting.\n");
	  return;
  }

  for(j = 0; j < TEST_SAMPLESIZE ; j++)
  {
    int szbl = sizeof(float);
    int sznxyz = sizeof(mynewxyz);
    //bdmpx_write(NULL, 4,sizeof(float),&backlerp,numVerts*sizeof(vec4_t),outXyz,numVerts*sizeof(vec4_t),outNormal,numVerts*8,newXyz);
    int ercd = bdmpx_read(NULL, 4, &szbl, &backlerp, 0, 0, 0, 0, &sznxyz, &mynewxyz[0]);
    
    if(ercd != 4)
      break;
    assert(szbl == sizeof(float));
    assert(sznxyz < sizeof(mynewxyz));

    surf.numVerts = sznxyz / 8;

    if(backlerp != 0)
    {
      int szoxyz = sizeof(myoldxyz);
      //bdmpx_write(NULL, 1,numVerts*8,oldXyz);
      ercd = bdmpx_read(NULL, 1, &szoxyz, &myoldxyz[0]);
      assert(szoxyz < sizeof(myoldxyz));
      assert(ercd == 1);
    }
	else
	{
		backlerp_was_zero++;
	}

	memset(&tess, 0, sizeof(tess));
	tess.guard = 0x5a5aa5a5;

	for (tested = 0, tess.numVertexes = 0; data_fn[tested].fp != 0; tested++, tess.numVertexes += SHADER_MAX_VERTEXES)
	{
		lerpmesh_fn callme = (lerpmesh_fn)data_fn[tested].fp;
		const char* myinfo = data_fn[tested].fname;

		cpu[tested].reset();
		timer[tested].reset();
		for (int r = 0; r < REPETITIONS; r++)
			callme(&surf, backlerp);
		elapsed[tested] = timer[tested].accum();
		cputime[tested] = cpu[tested].accum();
	}
    
	if (tess.guard != 0x5a5aa5a5)
	{
		printf("seg:%d guard failed %x\n", j, tess.guard);
	}
    float *scoord = tess.xyz[0], *snorm = tess.normal[0];
	float *vcoord, *vnorm;
	for (i = 1; i < tested; i++)
	{
		vcoord = tess.xyz[i * SHADER_MAX_VERTEXES];
		vnorm = tess.normal[i * SHADER_MAX_VERTEXES];
		for (k = 0; k < 4 * /*SHADER_MAX_VERTEXES*/surf.numVerts && err < MAX_ERRORS; k++)
		{
			if ((scoord[k] != vcoord[k]) || (snorm[k] != vnorm[k]))
			{
				printf("id:%d seg:%d samp:%d/%d scalar:%f|%f vect:%f|%f\n", i, j, k, surf.numVerts * 4, scoord[k], snorm[k], vcoord[k], vnorm[k]);
				err++;
			}
		}
	}
  }

  printf("LerpMesh Test segments %d (%d backlerps)\n", j, j-backlerp_was_zero);
  for (k = 0; k < tested; k++)
  {
	  const char* myinfo = data_fn[k].fname;
	  printf("%d %s: %4.4f %4.4f %4.4f %4.4f\n", k, myinfo, elapsed[k], elapsed[k] / elapsed[0], cputime[k], cputime[k] / cputime[0]);
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