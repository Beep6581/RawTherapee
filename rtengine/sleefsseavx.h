#include <immintrin.h>
#include <stdint.h>

#ifdef __SSE2__
#define VECTLENDP 2
#define VECTLENSP 4

typedef __m128d vdouble;
typedef __m128i vint;

typedef __m128 vfloat;
typedef __m128i vint2;
typedef __m128i vmask;

static vdouble vloadu(double *p) { return _mm_loadu_pd(p); }
static void vstoreu(double *p, vdouble v) { _mm_storeu_pd(p, v); }

static vfloat vloaduf(float *p) { return _mm_loadu_ps(p); }
static void vstoreuf(float *p, vfloat v) { _mm_storeu_ps(p, v); }

static vint2 vloadui2(int32_t *p) { return (vint2)_mm_loadu_si128((__m128i *)p); }
static void vstoreui2(int32_t *p, vint2 v) { _mm_storeu_si128((__m128i *)p, (__m128i)v); }
#endif

#ifdef ENABLE_AVX
#define VECTLENDP 4
#define VECTLENSP 8

typedef __m256d vdouble;
typedef __m128i vint;


typedef __m256 vfloat;
typedef struct {
  vint x, y;
} vint2;

static vdouble vloadu(double *p) { return _mm256_loadu_pd(p); }
static void vstoreu(double *p, vdouble v) { return _mm256_storeu_pd(p, v); }

static vfloat vloaduf(float *p) { return _mm256_loadu_ps(p); }
static void vstoreuf(float *p, vfloat v) { return _mm256_storeu_ps(p, v); }

static vint2 vloadui2(int32_t *p) {
  vint2 r;
  r.x = _mm_loadu_si128((__m128i *) p     );
  r.y = _mm_loadu_si128((__m128i *)(p + 4));
  return r;
}

static void vstoreui2(int32_t *p, vint2 v) {
  _mm_storeu_si128((__m128i *) p     , v.x);
  _mm_storeu_si128((__m128i *)(p + 4), v.y);
}
#endif

typedef struct {
  vdouble x, y;
} vdouble2;

vdouble xldexp(vdouble x, vint q);
vint xilogb(vdouble d);

vdouble xsin(vdouble d);
vdouble xcos(vdouble d);
vdouble2 xsincos(vdouble d);
vdouble xtan(vdouble d);
vdouble xasin(vdouble s);
vdouble xacos(vdouble s);
vdouble xatan(vdouble s);
vdouble xatan2(vdouble y, vdouble x);
vdouble xlog(vdouble d);
vdouble xexp(vdouble d);
vdouble xpow(vdouble x, vdouble y);

vdouble xsinh(vdouble d);
vdouble xcosh(vdouble d);
vdouble xtanh(vdouble d);
vdouble xasinh(vdouble s);
vdouble xacosh(vdouble s);
vdouble xatanh(vdouble s);

vdouble xcbrt(vdouble d);

vdouble xexp2(vdouble a);
vdouble xexp10(vdouble a);
vdouble xexpm1(vdouble a);
vdouble xlog10(vdouble a);
vdouble xlog1p(vdouble a);

//

typedef struct {
  vfloat x, y;
} vfloat2;

vfloat xsinf(vfloat d);
vfloat xcosf(vfloat d);
vfloat2 xsincosf(vfloat d);
vfloat xtanf(vfloat d);
vfloat xasinf(vfloat s);
vfloat xacosf(vfloat s);
vfloat xatanf(vfloat s);
vfloat xatan2f(vfloat y, vfloat x);
vfloat xlogf(vfloat d);
vfloat xlogf0(vfloat d);
vfloat xexpf(vfloat d);
vfloat xcbrtf(vfloat s);
