////////////////////////////////////////////////////////////////
//
//  this code was taken from http://shibatch.sourceforge.net/
//  Many thanks to the author: Naoki Shibata
//
////////////////////////////////////////////////////////////////
#ifndef __AVX__
#error Please specify -mavx.
#endif

#ifdef __GNUC__
#define INLINE __attribute__((always_inline))
#else
#define INLINE inline
#endif

#include <immintrin.h>
#include <stdint.h>

typedef __m256d vdouble;
typedef __m128i vint;
typedef __m256i vmask;

typedef __m256 vfloat;
typedef struct {
    vint x, y;
} vint2;

//

static INLINE vint vrint_vi_vd(vdouble vd)
{
    return _mm256_cvtpd_epi32(vd);
}
static INLINE vint vtruncate_vi_vd(vdouble vd)
{
    return _mm256_cvttpd_epi32(vd);
}
static INLINE vdouble vcast_vd_vi(vint vi)
{
    return _mm256_cvtepi32_pd(vi);
}
static INLINE vdouble vcast_vd_d(double d)
{
    return _mm256_set_pd(d, d, d, d);
}
static INLINE vint vcast_vi_i(int i)
{
    return _mm_set_epi32(i, i, i, i);
}

static INLINE vmask vreinterpret_vm_vd(vdouble vd)
{
    return (__m256i)vd;
}
static INLINE vdouble vreinterpret_vd_vm(vmask vm)
{
    return (__m256d)vm;
}

static INLINE vmask vreinterpret_vm_vf(vfloat vf)
{
    return (__m256i)vf;
}
static INLINE vfloat vreinterpret_vf_vm(vmask vm)
{
    return (__m256)vm;
}

//

static INLINE vfloat vcast_vf_f(float f)
{
    return _mm256_set_ps(f, f, f, f, f, f, f, f);
}

static INLINE vfloat vaddf(vfloat x, vfloat y)
{
    return _mm256_add_ps(x, y);
}
static INLINE vfloat vsubf(vfloat x, vfloat y)
{
    return _mm256_sub_ps(x, y);
}
static INLINE vfloat vmulf(vfloat x, vfloat y)
{
    return _mm256_mul_ps(x, y);
}
static INLINE vfloat vdivf(vfloat x, vfloat y)
{
    return _mm256_div_ps(x, y);
}
static INLINE vfloat vrecf(vfloat x)
{
    return vdivf(vcast_vf_f(1.0f), x);
}
static INLINE vfloat vsqrtf(vfloat x)
{
    return _mm256_sqrt_ps(x);
}
static INLINE vfloat vmaxf(vfloat x, vfloat y)
{
    return _mm256_max_ps(x, y);
}
static INLINE vfloat vminf(vfloat x, vfloat y)
{
    return _mm256_min_ps(x, y);
}

//

static INLINE vdouble vadd(vdouble x, vdouble y)
{
    return _mm256_add_pd(x, y);
}
static INLINE vdouble vsub(vdouble x, vdouble y)
{
    return _mm256_sub_pd(x, y);
}
static INLINE vdouble vmul(vdouble x, vdouble y)
{
    return _mm256_mul_pd(x, y);
}
static INLINE vdouble vdiv(vdouble x, vdouble y)
{
    return _mm256_div_pd(x, y);
}
static INLINE vdouble vrec(vdouble x)
{
    return _mm256_div_pd(_mm256_set_pd(1, 1, 1, 1), x);
}
static INLINE vdouble vsqrt(vdouble x)
{
    return _mm256_sqrt_pd(x);
}
static INLINE vdouble vmla(vdouble x, vdouble y, vdouble z)
{
    return vadd(vmul(x, y), z);
}

static INLINE vdouble vmax(vdouble x, vdouble y)
{
    return _mm256_max_pd(x, y);
}
static INLINE vdouble vmin(vdouble x, vdouble y)
{
    return _mm256_min_pd(x, y);
}

static INLINE vdouble vabs(vdouble d)
{
    return (__m256d)_mm256_andnot_pd(_mm256_set_pd(-0.0, -0.0, -0.0, -0.0), d);
}
static INLINE vdouble vneg(vdouble d)
{
    return (__m256d)_mm256_xor_pd(_mm256_set_pd(-0.0, -0.0, -0.0, -0.0), d);
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
    return _mm_slli_epi32 (x, c);
}
static INLINE vint vsrli(vint x, int c)
{
    return _mm_srli_epi32 (x, c);
}
static INLINE vint vsrai(vint x, int c)
{
    return _mm_srai_epi32 (x, c);
}

//

static INLINE vmask vandm(vmask x, vmask y)
{
    return (vmask)_mm256_and_pd((__m256d)x, (__m256d)y);
}
static INLINE vmask vandnotm(vmask x, vmask y)
{
    return (vmask)_mm256_andnot_pd((__m256d)x, (__m256d)y);
}
static INLINE vmask vorm(vmask x, vmask y)
{
    return (vmask)_mm256_or_pd((__m256d)x, (__m256d)y);
}
static INLINE vmask vxorm(vmask x, vmask y)
{
    return (vmask)_mm256_xor_pd((__m256d)x, (__m256d)y);
}

static INLINE vmask vmask_eq(vdouble x, vdouble y)
{
    return (__m256i)_mm256_cmp_pd(x, y, _CMP_EQ_OQ);
}
static INLINE vmask vmask_neq(vdouble x, vdouble y)
{
    return (__m256i)_mm256_cmp_pd(x, y, _CMP_NEQ_OQ);
}
static INLINE vmask vmask_lt(vdouble x, vdouble y)
{
    return (__m256i)_mm256_cmp_pd(x, y, _CMP_LT_OQ);
}
static INLINE vmask vmask_le(vdouble x, vdouble y)
{
    return (__m256i)_mm256_cmp_pd(x, y, _CMP_LE_OQ);
}
static INLINE vmask vmask_gt(vdouble x, vdouble y)
{
    return (__m256i)_mm256_cmp_pd(x, y, _CMP_GT_OQ);
}
static INLINE vmask vmask_ge(vdouble x, vdouble y)
{
    return (__m256i)_mm256_cmp_pd(x, y, _CMP_GE_OQ);
}

static INLINE vmask vmaskf_eq(vfloat x, vfloat y)
{
    return (__m256i)_mm256_cmp_ps(x, y, _CMP_EQ_OQ);
}
static INLINE vmask vmaskf_neq(vfloat x, vfloat y)
{
    return (__m256i)_mm256_cmp_ps(x, y, _CMP_NEQ_OQ);
}
static INLINE vmask vmaskf_lt(vfloat x, vfloat y)
{
    return (__m256i)_mm256_cmp_ps(x, y, _CMP_LT_OQ);
}
static INLINE vmask vmaskf_le(vfloat x, vfloat y)
{
    return (__m256i)_mm256_cmp_ps(x, y, _CMP_LE_OQ);
}
static INLINE vmask vmaskf_gt(vfloat x, vfloat y)
{
    return (__m256i)_mm256_cmp_ps(x, y, _CMP_GT_OQ);
}
static INLINE vmask vmaskf_ge(vfloat x, vfloat y)
{
    return (__m256i)_mm256_cmp_ps(x, y, _CMP_GE_OQ);
}

static INLINE vmask vmaski_eq(vint x, vint y)
{
    __m256d r = _mm256_cvtepi32_pd(_mm_and_si128(_mm_cmpeq_epi32(x, y), _mm_set_epi32(1, 1, 1, 1)));
    return vmask_eq(r, _mm256_set_pd(1, 1, 1, 1));
}

static INLINE vdouble vsel(vmask mask, vdouble x, vdouble y)
{
    return (__m256d)vorm(vandm(mask, (__m256i)x), vandnotm(mask, (__m256i)y));
}

static INLINE vint vseli_lt(vdouble d0, vdouble d1, vint x, vint y)
{
    __m128i mask = _mm256_cvtpd_epi32(_mm256_and_pd(_mm256_cmp_pd(d0, d1, _CMP_LT_OQ), _mm256_set_pd(1.0, 1.0, 1.0, 1.0)));
    mask = _mm_cmpeq_epi32(mask, _mm_set_epi32(1, 1, 1, 1));
    return vori(vandi(mask, x), vandnoti(mask, y));
}

//

static INLINE vint2 vcast_vi2_vm(vmask vm)
{
    vint2 r;
    r.x = _mm256_castsi256_si128(vm);
    r.y = _mm256_extractf128_si256(vm, 1);
    return r;
}

static INLINE vmask vcast_vm_vi2(vint2 vi)
{
    vmask m = _mm256_castsi128_si256(vi.x);
    m = _mm256_insertf128_si256(m, vi.y, 1);
    return m;
}

static INLINE vint2 vrint_vi2_vf(vfloat vf)
{
    return vcast_vi2_vm((vmask)_mm256_cvtps_epi32(vf));
}
static INLINE vint2 vtruncate_vi2_vf(vfloat vf)
{
    return vcast_vi2_vm((vmask)_mm256_cvttps_epi32(vf));
}
static INLINE vfloat vcast_vf_vi2(vint2 vi)
{
    return _mm256_cvtepi32_ps((vmask)vcast_vm_vi2(vi));
}
static INLINE vint2 vcast_vi2_i(int i)
{
    vint2 r;
    r.x = r.y = vcast_vi_i(i);
    return r;
}

static INLINE vint2 vaddi2(vint2 x, vint2 y)
{
    vint2 r;
    r.x = vaddi(x.x, y.x);
    r.y = vaddi(x.y, y.y);
    return r;
}
static INLINE vint2 vsubi2(vint2 x, vint2 y)
{
    vint2 r;
    r.x = vsubi(x.x, y.x);
    r.y = vsubi(x.y, y.y);
    return r;
}

static INLINE vint2 vandi2(vint2 x, vint2 y)
{
    vint2 r;
    r.x = vandi(x.x, y.x);
    r.y = vandi(x.y, y.y);
    return r;
}
static INLINE vint2 vandnoti2(vint2 x, vint2 y)
{
    vint2 r;
    r.x = vandnoti(x.x, y.x);
    r.y = vandnoti(x.y, y.y);
    return r;
}
static INLINE vint2 vori2(vint2 x, vint2 y)
{
    vint2 r;
    r.x = vori(x.x, y.x);
    r.y = vori(x.y, y.y);
    return r;
}
static INLINE vint2 vxori2(vint2 x, vint2 y)
{
    vint2 r;
    r.x = vxori(x.x, y.x);
    r.y = vxori(x.y, y.y);
    return r;
}

static INLINE vint2 vslli2(vint2 x, int c)
{
    vint2 r;
    r.x = vslli(x.x, c);
    r.y = vslli(x.y, c);
    return r;
}
static INLINE vint2 vsrli2(vint2 x, int c)
{
    vint2 r;
    r.x = vsrli(x.x, c);
    r.y = vsrli(x.y, c);
    return r;
}
static INLINE vint2 vsrai2(vint2 x, int c)
{
    vint2 r;
    r.x = vsrai(x.x, c);
    r.y = vsrai(x.y, c);
    return r;
}

static INLINE vmask vmaski2_eq(vint2 x, vint2 y)
{
    vint2 r;
    r.x = _mm_cmpeq_epi32(x.x, y.x);
    r.y = _mm_cmpeq_epi32(x.y, y.y);
    return vcast_vm_vi2(r);
}

static INLINE vint2 vseli2(vmask m, vint2 x, vint2 y)
{
    vint2 r, m2 = vcast_vi2_vm(m);
    r.x = vori(vandi(m2.x, x.x), vandnoti(m2.x, y.x));
    r.y = vori(vandi(m2.y, x.y), vandnoti(m2.y, y.y));
    return r;
}

//

static INLINE double vcast_d_vd(vdouble v)
{
    double s[4];
    _mm256_storeu_pd(s, v);
    return s[0];
}

static INLINE float vcast_f_vf(vfloat v)
{
    float s[8];
    _mm256_storeu_ps(s, v);
    return s[0];
}

static INLINE vmask vsignbit(vdouble d)
{
    return (vmask)_mm256_and_pd(d, _mm256_set_pd(-0.0, -0.0, -0.0, -0.0));
}

static INLINE vdouble vsign(vdouble d)
{
    return _mm256_or_pd(_mm256_set_pd(1.0, 1.0, 1.0, 1.0), (vdouble)vsignbit(d));
}

static INLINE vdouble vmulsign(vdouble x, vdouble y)
{
    return (__m256d)vxorm((__m256i)x, vsignbit(y));
}

static INLINE vmask vmask_isinf(vdouble d)
{
    return (vmask)_mm256_cmp_pd(vabs(d), _mm256_set_pd(INFINITY, INFINITY, INFINITY, INFINITY), _CMP_EQ_OQ);
}

static INLINE vmask vmask_ispinf(vdouble d)
{
    return (vmask)_mm256_cmp_pd(d, _mm256_set_pd(INFINITY, INFINITY, INFINITY, INFINITY), _CMP_EQ_OQ);
}

static INLINE vmask vmask_isminf(vdouble d)
{
    return (vmask)_mm256_cmp_pd(d, _mm256_set_pd(-INFINITY, -INFINITY, -INFINITY, -INFINITY), _CMP_EQ_OQ);
}

static INLINE vmask vmask_isnan(vdouble d)
{
    return (vmask)_mm256_cmp_pd(d, d, _CMP_NEQ_UQ);
}

static INLINE vdouble visinf(vdouble d)
{
    return _mm256_and_pd((vdouble)vmask_isinf(d), vsign(d));
}

static INLINE vdouble visinf2(vdouble d, vdouble m)
{
    return _mm256_and_pd((vdouble)vmask_isinf(d), _mm256_or_pd((vdouble)vsignbit(d), m));
}

static INLINE vdouble vpow2i(vint q)
{
    vint r;
    vdouble y;
    q = _mm_add_epi32(_mm_set_epi32(0x3ff, 0x3ff, 0x3ff, 0x3ff), q);
    q = _mm_slli_epi32(q, 20);
    r = (__m128i)_mm_shuffle_ps((__m128)q, (__m128)q, _MM_SHUFFLE(1, 0, 0, 0));
    y = _mm256_castpd128_pd256((__m128d)r);
    r = (__m128i)_mm_shuffle_ps((__m128)q, (__m128)q, _MM_SHUFFLE(3, 2, 2, 2));
    y = _mm256_insertf128_pd(y, (__m128d)r, 1);
    y = _mm256_and_pd(y, (__m256d)_mm256_set_epi32(0xfff00000, 0, 0xfff00000, 0, 0xfff00000, 0, 0xfff00000, 0));
    return y;
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
    vint q, r, c;
    vmask m = vmask_lt(d, vcast_vd_d(4.9090934652977266E-91));
    d = vsel(m, vmul(vcast_vd_d(2.037035976334486E90), d), d);
    c = _mm256_cvtpd_epi32(vsel(m, vcast_vd_d(300 + 0x3fe), vcast_vd_d(0x3fe)));
    q = (__m128i)_mm256_castpd256_pd128(d);
    q = (__m128i)_mm_shuffle_ps((__m128)q, _mm_set_ps(0, 0, 0, 0), _MM_SHUFFLE(0, 0, 3, 1));
    r = (__m128i)_mm256_extractf128_pd(d, 1);
    r = (__m128i)_mm_shuffle_ps(_mm_set_ps(0, 0, 0, 0), (__m128)r, _MM_SHUFFLE(3, 1, 0, 0));
    q = _mm_or_si128(q, r);
    q = _mm_srli_epi32(q, 20);
    q = _mm_sub_epi32(q, c);
    return q;
}

static INLINE vdouble vupper(vdouble d)
{
    return (__m256d)_mm256_and_pd(d, (vdouble)_mm256_set_epi32(0xffffffff, 0xf8000000, 0xffffffff, 0xf8000000, 0xffffffff, 0xf8000000, 0xffffffff, 0xf8000000));
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
    return dd((__m256d)vorm(vandm(mask, (__m256i)x.x), vandnotm(mask, (__m256i)y.x)),
              (__m256d)vorm(vandm(mask, (__m256i)x.y), vandnotm(mask, (__m256i)y.y)));
}

static INLINE vdouble2 abs_d(vdouble2 x)
{
    return dd((__m256d)_mm256_xor_pd(_mm256_and_pd(_mm256_set_pd(-0.0, -0.0, -0.0, -0.0), x.x), x.x),
              (__m256d)_mm256_xor_pd(_mm256_and_pd(_mm256_set_pd(-0.0, -0.0, -0.0, -0.0), x.x), x.y));
}
