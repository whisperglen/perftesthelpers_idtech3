
#ifndef PLATFORM_H
#define PLATFORM_H

#define QDECL

#define id386 0
#define idx64 0
#define idsse 0
#define idarm 0
#define idneon 0

#ifdef _WIN32

#undef QDECL
#define QDECL __cdecl

#define ID_INLINE __inline
#define QALIGNA(x) __declspec(align(x))
#define QALIGNB(x)

#if defined( _M_IX86 )
#if _M_IX86_FP >= 2
#undef idsse
#define idsse 1
#endif
#undef id386
#define id386 1
#endif

#if defined( _M_AMD64 )
#undef idx64
#define idx64 1
//#undef id386
//#define id386 1
#undef idsse
#define idsse 1
#endif

#endif

#ifdef __linux__

#define ID_INLINE inline
#define QALIGNA(x)
#define QALIGNB(x) __attribute__((aligned(x)))

#if defined __i386__
#if defined __SSE__
#undef idsse
#define idsse 1
#endif
#undef id386
#define id386 1
#endif

#if defined __x86_64__
#undef idx64
#define idx64 1
#undef idsse
#define idsse 1
#endif

#if defined __arm__
#if defined __ARM_NEON__
#undef idneon
#define idneon 1
#endif
#undef idarm
#define idarm 1
#endif

#endif

#endif