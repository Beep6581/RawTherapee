#ifndef __SSE2__
#error Please specify -msse2.
#endif

#ifdef __GNUC__
#define INLINE __inline
//#define INLINE __attribute__((always_inline))
#else
#define INLINE inline
#endif

#include <x86intrin.h>

#include <stdint.h>

typedef __m128d vdouble;
typedef __m128i vint;
typedef __m128i vmask;

typedef __m128 vfloat;
typedef __m128i vint2;

//
#ifdef __GNUC__
#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 9) || __GNUC__ > 4
#define LVF(x) _mm_load_ps(&x)
#define LVFU(x) _mm_loadu_ps(&x)
#define STVF(x,y) _mm_store_ps(&x,y)
#define STVFU(x,y) _mm_storeu_ps(&x,y)
#else // there is a bug in gcc 4.7.x when using openmp and aligned memory and -O3
#define LVF(x) _mm_loadu_ps(&x)
#define LVFU(x) _mm_loadu_ps(&x)
#define STVF(x,y) _mm_storeu_ps(&x,y)
#define STVFU(x,y) _mm_storeu_ps(&x,y)
#endif
#else
#define LVF(x) _mm_load_ps(&x)
#define LVFU(x) _mm_loadu_ps(&x)
#define STVF(x,y) _mm_store_ps(&x,y)
#define STVFU(x,y) _mm_storeu_ps(&x,y)
#endif

// Load 8 floats from a and combine a[0],a[2],a[4] and a[6] into a vector of 4 floats
#define LC2VFU(a) _mm_shuffle_ps( LVFU(a),  _mm_loadu_ps(  (&a) + 4 ), _MM_SHUFFLE( 2,0,2,0 ) )

// Store a vector of 4 floats in a[0],a[2],a[4] and a[6]
#if defined(__x86_64__) && defined(__SSE4_1__)
// SSE4.1 => use _mm_blend_ps instead of _mm_set_epi32 and vself
#define STC2VFU(a,v) {\
                         __m128 TST1V = _mm_loadu_ps(&a);\
                         __m128 TST2V = _mm_shuffle_ps(v,v,_MM_SHUFFLE( 1,1,0,0 ));\
                         _mm_storeu_ps(&a, _mm_blend_ps(TST1V,TST2V,5));\
                         TST1V = _mm_loadu_ps((&a)+4);\
                         TST2V = _mm_shuffle_ps(v,v,_MM_SHUFFLE( 3,3,2,2 ));\
                         _mm_storeu_ps((&a)+4, _mm_blend_ps(TST1V,TST2V,5));\
                     }
#else
#define STC2VFU(a,v) {\
                         __m128 TST1V = _mm_loadu_ps(&a);\
                         __m128 TST2V = _mm_shuffle_ps(v,v,_MM_SHUFFLE( 1,1,0,0 ));\
                         vmask cmask = _mm_set_epi32(0xffffffff,0,0xffffffff,0);\
                         _mm_storeu_ps(&a, vself(cmask,TST1V,TST2V));\
                         TST1V = _mm_loadu_ps((&a)+4);\
                         TST2V = _mm_shuffle_ps(v,v,_MM_SHUFFLE( 3,3,2,2 ));\
                         _mm_storeu_ps((&a)+4, vself(cmask,TST1V,TST2V));\
                     }
#endif

#define ZEROV _mm_setzero_ps()
#define F2V(a) _mm_set1_ps((a))

static INLINE vint vrint_vi_vd(vdouble vd)
{
    return _mm_cvtpd_epi32(vd);
}
static INLINE vint vtruncate_vi_vd(vdouble vd)
{
    return _mm_cvttpd_epi32(vd);
}
static INLINE vdouble vcast_vd_vi(vint vi)
{
    return _mm_cvtepi32_pd(vi);
}
static INLINE vdouble vcast_vd_d(double d)
{
    return _mm_set_pd(d, d);
}
static INLINE vint vcast_vi_i(int i)
{
    return _mm_set_epi32(0, 0, i, i);
}

static INLINE vmask vreinterpret_vm_vd(vdouble vd)
{
    return (__m128i)vd;
}
static INLINE vdouble vreinterpret_vd_vm(vint vm)
{
    return (__m128d)vm;
}

static INLINE vmask vreinterpret_vm_vf(vfloat vf)
{
    return (__m128i)vf;
}
static INLINE vfloat vreinterpret_vf_vm(vmask vm)
{
    return (__m128)vm;
}

//

static INLINE vfloat vcast_vf_f(float f)
{
    return _mm_set_ps(f, f, f, f);
}

// Don't use intrinsics here. Newer gcc versions (>= 4.9, maybe also before 4.9) generate better code when not using intrinsics
// example: vaddf(vmulf(a,b),c) will generate an FMA instruction when build for chips with that feature only when vaddf and vmulf don't use intrinsics
static INLINE vfloat vaddf(vfloat x, vfloat y)
{
    return x + y;
}
static INLINE vfloat vsubf(vfloat x, vfloat y)
{
    return x - y;
}
static INLINE vfloat vmulf(vfloat x, vfloat y)
{
    return x * y;
}
static INLINE vfloat vdivf(vfloat x, vfloat y)
{
    return x / y;
}
// Also don't use intrinsic here: Some chips support FMA instructions with 3 and 4 operands
// 3 operands: a = a*b+c, b = a*b+c, c = a*b+c // destination has to be one of a,b,c
// 4 operands: d = a*b+c // destination does not have to be one of a,b,c
// gcc will use the one which fits best when not using intrinsics. With using intrinsics that's not possible
static INLINE vfloat vmlaf(vfloat x, vfloat y, vfloat z) {
    return x * y + z;
}
static INLINE vfloat vrecf(vfloat x)
{
    return vdivf(vcast_vf_f(1.0f), x);
}
static INLINE vfloat vsqrtf(vfloat x)
{
    return _mm_sqrt_ps(x);
}
static INLINE vfloat vmaxf(vfloat x, vfloat y)
{
    return _mm_max_ps(x, y);
}
static INLINE vfloat vminf(vfloat x, vfloat y)
{
    return _mm_min_ps(x, y);
}

//

static INLINE vdouble vadd(vdouble x, vdouble y)
{
    return _mm_add_pd(x, y);
}
static INLINE vdouble vsub(vdouble x, vdouble y)
{
    return _mm_sub_pd(x, y);
}
static INLINE vdouble vmul(vdouble x, vdouble y)
{
    return _mm_mul_pd(x, y);
}
static INLINE vdouble vdiv(vdouble x, vdouble y)
{
    return _mm_div_pd(x, y);
}
static INLINE vdouble vrec(vdouble x)
{
    return _mm_div_pd(_mm_set_pd(1, 1), x);
}
static INLINE vdouble vsqrt(vdouble x)
{
    return _mm_sqrt_pd(x);
}
static INLINE vdouble vmla(vdouble x, vdouble y, vdouble z)
{
    return vadd(vmul(x, y), z);
}

static INLINE vdouble vmax(vdouble x, vdouble y)
{
    return _mm_max_pd(x, y);
}
static INLINE vdouble vmin(vdouble x, vdouble y)
{
    return _mm_min_pd(x, y);
}

static INLINE vdouble vabs(vdouble d)
{
    return (__m128d)_mm_andnot_pd(_mm_set_pd(-0.0, -0.0), d);
}
static INLINE vdouble vneg(vdouble d)
{
    return (__m128d)_mm_xor_pd(_mm_set_pd(-0.0, -0.0), d);
}

//

static INLINE vint vaddi(vint x, vint y)
{
    return _mm_add_epi32(x, y);
}
static INLINE vint vsubi(vint x, vint y)
{
    return _mm_sub_epi32(x, y);
}

static INLINE vint vandi(vint x, vint y)
{
    return _mm_and_si128(x, y);
}
static INLINE vint vandnoti(vint x, vint y)
{
    return _mm_andnot_si128(x, y);
}
static INLINE vint vori(vint x, vint y)
{
    return _mm_or_si128(x, y);
}
static INLINE vint vxori(vint x, vint y)
{
    return _mm_xor_si128(x, y);
}

static INLINE vint vslli(vint x, int c)
{
    return _mm_slli_epi32(x, c);
}
static INLINE vint vsrli(vint x, int c)
{
    return _mm_srli_epi32(x, c);
}
static INLINE vint vsrai(vint x, int c)
{
    return _mm_srai_epi32(x, c);
}

//

static INLINE vmask vandm(vmask x, vmask y)
{
    return _mm_and_si128(x, y);
}
static INLINE vmask vandnotm(vmask x, vmask y)
{
    return _mm_andnot_si128(x, y);
}
static INLINE vmask vorm(vmask x, vmask y)
{
    return _mm_or_si128(x, y);
}
static INLINE vmask vxorm(vmask x, vmask y)
{
    return _mm_xor_si128(x, y);
}
static INLINE vmask vnotm(vmask x)
{
    return _mm_xor_si128(x, _mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128()));
}

static INLINE vmask vmask_eq(vdouble x, vdouble y)
{
    return (__m128i)_mm_cmpeq_pd(x, y);
}
static INLINE vmask vmask_neq(vdouble x, vdouble y)
{
    return (__m128i)_mm_cmpneq_pd(x, y);
}
static INLINE vmask vmask_lt(vdouble x, vdouble y)
{
    return (__m128i)_mm_cmplt_pd(x, y);
}
static INLINE vmask vmask_le(vdouble x, vdouble y)
{
    return (__m128i)_mm_cmple_pd(x, y);
}
static INLINE vmask vmask_gt(vdouble x, vdouble y)
{
    return (__m128i)_mm_cmpgt_pd(x, y);
}
static INLINE vmask vmask_ge(vdouble x, vdouble y)
{
    return (__m128i)_mm_cmpge_pd(x, y);
}

static INLINE vmask vmaskf_eq(vfloat x, vfloat y)
{
    return (__m128i)_mm_cmpeq_ps(x, y);
}
static INLINE vmask vmaskf_neq(vfloat x, vfloat y)
{
    return (__m128i)_mm_cmpneq_ps(x, y);
}
static INLINE vmask vmaskf_lt(vfloat x, vfloat y)
{
    return (__m128i)_mm_cmplt_ps(x, y);
}
static INLINE vmask vmaskf_le(vfloat x, vfloat y)
{
    return (__m128i)_mm_cmple_ps(x, y);
}
static INLINE vmask vmaskf_gt(vfloat x, vfloat y)
{
    return (__m128i)_mm_cmpgt_ps(x, y);
}
static INLINE vmask vmaskf_ge(vfloat x, vfloat y)
{
    return (__m128i)_mm_cmpge_ps(x, y);
}


static INLINE vmask vmaski_eq(vint x, vint y)
{
    __m128 s = (__m128)_mm_cmpeq_epi32(x, y);
    return (__m128i)_mm_shuffle_ps(s, s, _MM_SHUFFLE(1, 1, 0, 0));
}

static INLINE vdouble vsel(vmask mask, vdouble x, vdouble y)
{
    return (__m128d)vorm(vandm(mask, (__m128i)x), vandnotm(mask, (__m128i)y));
}

static INLINE vint vseli_lt(vdouble d0, vdouble d1, vint x, vint y)
{
    vmask mask = (vmask)_mm_cmpeq_ps(_mm_cvtpd_ps((vdouble)vmask_lt(d0, d1)), _mm_set_ps(0, 0, 0, 0));
    return vori(vandnoti(mask, x), vandi(mask, y));
}

//

static INLINE vint2 vcast_vi2_vm(vmask vm)
{
    return (vint2)vm;
}
static INLINE vmask vcast_vm_vi2(vint2 vi)
{
    return (vmask)vi;
}

static INLINE vint2 vrint_vi2_vf(vfloat vf)
{
    return _mm_cvtps_epi32(vf);
}
static INLINE vint2 vtruncate_vi2_vf(vfloat vf)
{
    return _mm_cvttps_epi32(vf);
}
static INLINE vfloat vcast_vf_vi2(vint2 vi)
{
    return _mm_cvtepi32_ps(vcast_vm_vi2(vi));
}
static INLINE vint2 vcast_vi2_i(int i)
{
    return _mm_set_epi32(i, i, i, i);
}

static INLINE vint2 vaddi2(vint2 x, vint2 y)
{
    return vaddi(x, y);
}
static INLINE vint2 vsubi2(vint2 x, vint2 y)
{
    return vsubi(x, y);
}

static INLINE vint2 vandi2(vint2 x, vint2 y)
{
    return vandi(x, y);
}
static INLINE vint2 vandnoti2(vint2 x, vint2 y)
{
    return vandnoti(x, y);
}
static INLINE vint2 vori2(vint2 x, vint2 y)
{
    return vori(x, y);
}
static INLINE vint2 vxori2(vint2 x, vint2 y)
{
    return vxori(x, y);
}

static INLINE vint2 vslli2(vint2 x, int c)
{
    return vslli(x, c);
}
static INLINE vint2 vsrli2(vint2 x, int c)
{
    return vsrli(x, c);
}
static INLINE vint2 vsrai2(vint2 x, int c)
{
    return vsrai(x, c);
}

static INLINE vmask vmaski2_eq(vint2 x, vint2 y)
{
    return _mm_cmpeq_epi32(x, y);
}
static INLINE vint2 vseli2(vmask m, vint2 x, vint2 y)
{
    return vorm(vandm(m, x), vandnotm(m, y));
}

//

static INLINE double vcast_d_vd(vdouble v)
{
    double s[2];
    _mm_storeu_pd(s, v);
    return s[0];
}

static INLINE float vcast_f_vf(vfloat v)
{
    float s[4];
    _mm_storeu_ps(s, v);
    return s[0];
}

static INLINE vmask vsignbit(vdouble d)
{
    return _mm_and_si128((__m128i)d, _mm_set_epi32(0x80000000, 0x0, 0x80000000, 0x0));
}

static INLINE vdouble vsign(vdouble d)
{
    return (__m128d)_mm_or_si128((__m128i)_mm_set_pd(1, 1), _mm_and_si128((__m128i)d, _mm_set_epi32(0x80000000, 0x0, 0x80000000, 0x0)));
}

static INLINE vdouble vmulsign(vdouble x, vdouble y)
{
    return (__m128d)vxori((__m128i)x, vsignbit(y));
}

static INLINE vmask vmask_isinf(vdouble d)
{
    return (vmask)_mm_cmpeq_pd(vabs(d), _mm_set_pd(INFINITY, INFINITY));
}

static INLINE vmask vmask_ispinf(vdouble d)
{
    return (vmask)_mm_cmpeq_pd(d, _mm_set_pd(INFINITY, INFINITY));
}

static INLINE vmask vmask_isminf(vdouble d)
{
    return (vmask)_mm_cmpeq_pd(d, _mm_set_pd(-INFINITY, -INFINITY));
}

static INLINE vmask vmask_isnan(vdouble d)
{
    return (vmask)_mm_cmpneq_pd(d, d);
}

static INLINE vdouble visinf(vdouble d)
{
    return (__m128d)_mm_and_si128(vmask_isinf(d), _mm_or_si128(vsignbit(d), (__m128i)_mm_set_pd(1, 1)));
}

static INLINE vdouble visinf2(vdouble d, vdouble m)
{
    return (__m128d)_mm_and_si128(vmask_isinf(d), _mm_or_si128(vsignbit(d), (__m128i)m));
}

//

static INLINE vdouble vpow2i(vint q)
{
    q = _mm_add_epi32(_mm_set_epi32(0x0, 0x0, 0x3ff, 0x3ff), q);
    q = (__m128i)_mm_shuffle_ps((__m128)q, (__m128)q, _MM_SHUFFLE(1, 3, 0, 3));
    return (__m128d)_mm_slli_epi32(q, 20);
}

static INLINE vdouble vldexp(vdouble x, vint q)
{
    vint m = _mm_srai_epi32(q, 31);
    m = _mm_slli_epi32(_mm_sub_epi32(_mm_srai_epi32(_mm_add_epi32(m, q), 9), m), 7);
    q = _mm_sub_epi32(q, _mm_slli_epi32(m, 2));
    vdouble y = vpow2i(m);
    return vmul(vmul(vmul(vmul(vmul(x, y), y), y), y), vpow2i(q));
}

static INLINE vint vilogbp1(vdouble d)
{
    vint m = vmask_lt(d, vcast_vd_d(4.9090934652977266E-91));
    d = vsel(m, vmul(vcast_vd_d(2.037035976334486E90), d), d);
    __m128i q = _mm_and_si128((__m128i)d, _mm_set_epi32(((1 << 12) - 1) << 20, 0, ((1 << 12) - 1) << 20, 0));
    q = _mm_srli_epi32(q, 20);
    q = vorm(vandm   (m, _mm_sub_epi32(q, _mm_set_epi32(300 + 0x3fe, 0, 300 + 0x3fe, 0))),
             vandnotm(m, _mm_sub_epi32(q, _mm_set_epi32(      0x3fe, 0,       0x3fe, 0))));
    q = (__m128i)_mm_shuffle_ps((__m128)q, (__m128)q, _MM_SHUFFLE(0, 0, 3, 1));
    return q;
}

static INLINE vdouble vupper(vdouble d)
{
    return (__m128d)_mm_and_si128((__m128i)d, _mm_set_epi32(0xffffffff, 0xf8000000, 0xffffffff, 0xf8000000));
}

//

typedef struct {
    vdouble x, y;
} vdouble2;

static INLINE vdouble2 dd(vdouble h, vdouble l)
{
    vdouble2 ret = {h, l};
    return ret;
}

static INLINE vdouble2 vsel2(vmask mask, vdouble2 x, vdouble2 y)
{
    return dd((__m128d)vorm(vandm(mask, (__m128i)x.x), vandnotm(mask, (__m128i)y.x)),
              (__m128d)vorm(vandm(mask, (__m128i)x.y), vandnotm(mask, (__m128i)y.y)));
}

static INLINE vdouble2 abs_d(vdouble2 x)
{
    return dd((__m128d)_mm_xor_pd(_mm_and_pd(_mm_set_pd(-0.0, -0.0), x.x), x.x),
              (__m128d)_mm_xor_pd(_mm_and_pd(_mm_set_pd(-0.0, -0.0), x.x), x.y));
}
