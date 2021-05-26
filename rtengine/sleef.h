////////////////////////////////////////////////////////////////
//
//  this code was taken from http://shibatch.sourceforge.net/
//  Many thanks to the author of original version: Naoki Shibata
//
//   Copyright Naoki Shibata and contributors 2010 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file sleef_LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//
//  This version contains modifications made by Ingo Weyrich
//
////////////////////////////////////////////////////////////////
#pragma once

#include <assert.h>
#include <stdint.h>
#include "rt_math.h"
#include "opthelper.h"

//
// Double precision functions (see sleef-2.121/purec/sleefdp.c)
// All functions available, not all used
//

#define PI_A 3.1415926218032836914
#define PI_B 3.1786509424591713469e-08
#define PI_C 1.2246467864107188502e-16
#define PI_D 1.2736634327021899816e-24

#define M_2_PI_H 0.63661977236758138243
#define M_2_PI_L -3.9357353350364971764e-17

#define TRIGRANGEMAX 1e+14
#define SQRT_DBL_MAX 1.3407807929942596355e+154

#define L2U .69314718055966295651160180568695068359375
#define L2L .28235290563031577122588448175013436025525412068e-12
#define R_LN2 1.442695040888963407359924681001892137426645954152985934135449406931

__inline int64_t doubleToRawLongBits(double d) {
    union {
        double f;
        int64_t i;
    } tmp;
    tmp.f = d;
    return tmp.i;
}

__inline double longBitsToDouble(int64_t i) {
    union {
        double f;
        int64_t i;
    } tmp;
    tmp.i = i;
    return tmp.f;
}

__inline double xfabs(double x) {
    return longBitsToDouble(0x7fffffffffffffffLL & doubleToRawLongBits(x));
}

__inline double mulsign(double x, double y) {
    return longBitsToDouble(doubleToRawLongBits(x) ^ (doubleToRawLongBits(y) & (1LL << 63)));
}

__inline double sign(double d) { return mulsign(1, d); }
__inline double mla(double x, double y, double z) { return x * y + z; }
__inline double xrint(double x) { return x < 0 ? (int)(x - 0.5) : (int)(x + 0.5); }
__inline double xtrunc(double x) { return (double)(int)x; }

__inline int xisnan(double x) { return x != x; }
__inline int xisinf(double x) { return x == rtengine::RT_INFINITY || x == -rtengine::RT_INFINITY; }
__inline int xisminf(double x) { return x == -rtengine::RT_INFINITY; }
__inline int xispinf(double x) { return x == rtengine::RT_INFINITY; }
__inline int xisnegzero(double x) { return doubleToRawLongBits(x) == doubleToRawLongBits(-0.0); }

__inline int xisint(double d) {
  double x = d - (double)(1 << 31) * (int)(d * (1.0 / (1 << 31)));
  return (x == (int)x) || (xfabs(d) >= (double)(1LL << 52));
}

__inline double pow2i(int q) {
  return longBitsToDouble(((int64_t)(q + 0x3ff)) << 52);
}

__inline double ldexpk(double x, int q) {
    double u;
    int m;
    m = q >> 31;
    m = (((m + q) >> 9) - m) << 7;
    q = q - (m << 2);
    m += 0x3ff;
    m = m < 0     ? 0     : m;
    m = m > 0x7ff ? 0x7ff : m;
    u = longBitsToDouble(((int64_t)m) << 52);
    x = x * u * u * u * u;
    u = longBitsToDouble(((int64_t)(q + 0x3ff)) << 52);
    return x * u;
}

__inline double xldexp(double x, int q) { return ldexpk(x, q); }

__inline int ilogbk(double d) {
    int m = d < 4.9090934652977266E-91;
    d = m ? 2.037035976334486E90 * d : d;
    int q = (doubleToRawLongBits(d) >> 52) & 0x7ff;
    q = m ? q - (300 + 0x03ff) : q - 0x03ff;
    return q;
}

__inline int xilogb(double d) {
    int e = ilogbk(xfabs(d));
    e = d == 0.0  ? FP_ILOGB0 : e;
    e = xisnan(d) ? FP_ILOGBNAN : e;
    e = xisinf(d) ? INT_MAX : e;
    return e;
}

typedef struct {
    double x, y;
} double2;

typedef struct {
    float x, y;
} float2;

__inline double upper(double d) {
    return longBitsToDouble(doubleToRawLongBits(d) & 0xfffffffff8000000LL);
}

__inline double2 dd(double h, double l) {
    double2 ret;
    ret.x = h; ret.y = l;
    return ret;
}

__inline double2 ddnormalize_d2_d2(double2 t) {
    double2 s;

    s.x = t.x + t.y;
    s.y = t.x - s.x + t.y;

    return s;
}

__inline double2 ddscale_d2_d2_d(double2 d, double s) {
    double2 r;

    r.x = d.x * s;
    r.y = d.y * s;

    return r;
}

__inline double2 ddneg_d2_d2(double2 d) {
    double2 r;
    
    r.x = -d.x;
    r.y = -d.y;
    
    return r;
}

__inline double2 ddadd_d2_d_d(double x, double y) {
    // |x| >= |y|
    
    double2 r;
    
    // Snip assertion
    
    r.x = x + y;
    r.y = x - r.x + y;
    
    return r;
}

__inline double2 ddadd2_d2_d_d(double x, double y) {
    double2 r;

    r.x = x + y;
    double v = r.x - x;
    r.y = (x - (r.x - v)) + (y - v);

    return r;
}

__inline double2 ddadd_d2_d2_d(double2 x, double y) {
    // |x| >= |y|

    double2 r;

    // Snip assertion

    r.x = x.x + y;
    r.y = x.x - r.x + y + x.y;

    return r;
}

__inline double2 ddadd2_d2_d2_d(double2 x, double y) {
    // |x| >= |y|

    double2 r;

    r.x  = x.x + y;
    double v = r.x - x.x;
    r.y = (x.x - (r.x - v)) + (y - v);
    r.y += x.y;

    return r;
}

__inline double2 ddadd_d2_d_d2(double x, double2 y) {
    // |x| >= |y|

    double2 r;

    // Snip assertion

    r.x = x + y.x;
    r.y = x - r.x + y.x + y.y;

    return r;
}

__inline double2 ddadd2_d2_d_d2(double x, double2 y) {
    double2 r;
    
    r.x  = x + y.x;
    double v = r.x - x;
    r.y = (x - (r.x - v)) + (y.x - v) + y.y;
    
    return r;
}

__inline double2 ddadd_d2_d2_d2(double2 x, double2 y) {
    // |x| >= |y|

    double2 r;

    // Snip assertion

    r.x = x.x + y.x;
    r.y = x.x - r.x + y.x + x.y + y.y;

    return r;
}

__inline double2 ddadd2_d2_d2_d2(double2 x, double2 y) {
    double2 r;

    r.x  = x.x + y.x;
    double v = r.x - x.x;
    r.y = (x.x - (r.x - v)) + (y.x - v);
    r.y += x.y + y.y;

    return r;
}

__inline double2 ddsub_d2_d2_d2(double2 x, double2 y) {
    // |x| >= |y|
    
    double2 r;
    
    // Snip assertion
    
    r.x = x.x - y.x;
    r.y = x.x - r.x - y.x + x.y - y.y;
    
    return r;
}

__inline double2 dddiv_d2_d2_d2(double2 n, double2 d) {
    double t = 1.0 / d.x;
    double dh  = upper(d.x), dl  = d.x - dh;
    double th  = upper(t  ), tl  = t   - th;
    double nhh = upper(n.x), nhl = n.x - nhh;

    double2 q;

    q.x = n.x * t;

    double u = -q.x + nhh * th + nhh * tl + nhl * th + nhl * tl +
        q.x * (1 - dh * th - dh * tl - dl * th - dl * tl);

    q.y = t * (n.y - q.x * d.y) + u;

    return q;
}

__inline double2 ddmul_d2_d_d(double x, double y) {
    double xh = upper(x), xl = x - xh;
    double yh = upper(y), yl = y - yh;
    double2 r;

    r.x = x * y;
    r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl;

    return r;
}

__inline double2 ddmul_d2_d2_d(double2 x, double y) {
    double xh = upper(x.x), xl = x.x - xh;
    double yh = upper(y  ), yl = y   - yh;
    double2 r;

    r.x = x.x * y;
    r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl + x.y * y;

    return r;
}

__inline double2 ddmul_d2_d2_d2(double2 x, double2 y) {
    double xh = upper(x.x), xl = x.x - xh;
    double yh = upper(y.x), yl = y.x - yh;
    double2 r;

    r.x = x.x * y.x;
    r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl + x.x * y.y + x.y * y.x;

    return r;
}

__inline double2 ddsqu_d2_d2(double2 x) {
    double xh = upper(x.x), xl = x.x - xh;
    double2 r;

    r.x = x.x * x.x;
    r.y = xh * xh - r.x + (xh + xh) * xl + xl * xl + x.x * (x.y + x.y);

    return r;
}

__inline double2 ddrec_d2_d(double d) {
    double t = 1.0 / d;
    double dh = upper(d), dl = d - dh;
    double th = upper(t), tl = t - th;
    double2 q;

    q.x = t;
    q.y = t * (1 - dh * th - dh * tl - dl * th - dl * tl);

    return q;
}

__inline double2 ddrec_d2_d2(double2 d) {
    double t = 1.0 / d.x;
    double dh = upper(d.x), dl = d.x - dh;
    double th = upper(t  ), tl = t   - th;
    double2 q;
    
    q.x = t;
    q.y = t * (1 - dh * th - dh * tl - dl * th - dl * tl - d.y * t);
    
    return q;
}

__inline double2 ddsqrt_d2_d2(double2 d) {
    double t = sqrt(d.x + d.y);
    return ddscale_d2_d2_d(ddmul_d2_d2_d2(ddadd2_d2_d2_d2(d, ddmul_d2_d_d(t, t)), ddrec_d2_d(t)), 0.5);
}

__inline double atan2k(double y, double x) {
    double s, t, u;
    int q = 0;

    if (x < 0) { x = -x; q = -2; }
    if (y > x) { t = x; x = y; y = -t; q += 1; }

    s = y / x;
    t = s * s;

    u = -1.88796008463073496563746e-05;
    u = u * t + (0.000209850076645816976906797);
    u = u * t + (-0.00110611831486672482563471);
    u = u * t + (0.00370026744188713119232403);
    u = u * t + (-0.00889896195887655491740809);
    u = u * t + (0.016599329773529201970117);
    u = u * t + (-0.0254517624932312641616861);
    u = u * t + (0.0337852580001353069993897);
    u = u * t + (-0.0407629191276836500001934);
    u = u * t + (0.0466667150077840625632675);
    u = u * t + (-0.0523674852303482457616113);
    u = u * t + (0.0587666392926673580854313);
    u = u * t + (-0.0666573579361080525984562);
    u = u * t + (0.0769219538311769618355029);
    u = u * t + (-0.090908995008245008229153);
    u = u * t + (0.111111105648261418443745);
    u = u * t + (-0.14285714266771329383765);
    u = u * t + (0.199999999996591265594148);
    u = u * t + (-0.333333333333311110369124);

    t = u * t * s + s;
    t = q * (rtengine::RT_PI_2) + t;

    return t;
}

__inline double xatan2(double y, double x) {
    double r = atan2k(xfabs(y), x);

    r = mulsign(r, x);
    if (xisinf(x) || x == 0) r = rtengine::RT_PI_2 - (xisinf(x) ? (sign(x) * (rtengine::RT_PI_2)) : 0);
    if (xisinf(y)          ) r = rtengine::RT_PI_2 - (xisinf(x) ? (sign(x) * (rtengine::RT_PI*1/4)) : 0);
    if (             y == 0) r = (sign(x) == -1 ? rtengine::RT_PI : 0);

    return xisnan(x) || xisnan(y) ? rtengine::RT_NAN : mulsign(r, y);
}

__inline double xasin(double d) {
    return mulsign(atan2k(xfabs(d), sqrt((1+d)*(1-d))), d);
}

__inline double xacos(double d) {
    return mulsign(atan2k(sqrt((1+d)*(1-d)), xfabs(d)), d) + (sign(d) == -1 ? rtengine::RT_PI : 0);
}

__inline double xatan(double s) {
    double t, u;
    int q = 0;

    if (sign(s) == -1) { s = -s; q = 2; }
    if (s > 1) { s = 1.0 / s; q |= 1; }

    t = s * s;

    u = -1.88796008463073496563746e-05;
    u = u * t + (0.000209850076645816976906797);
    u = u * t + (-0.00110611831486672482563471);
    u = u * t + (0.00370026744188713119232403);
    u = u * t + (-0.00889896195887655491740809);
    u = u * t + (0.016599329773529201970117);
    u = u * t + (-0.0254517624932312641616861);
    u = u * t + (0.0337852580001353069993897);
    u = u * t + (-0.0407629191276836500001934);
    u = u * t + (0.0466667150077840625632675);
    u = u * t + (-0.0523674852303482457616113);
    u = u * t + (0.0587666392926673580854313);
    u = u * t + (-0.0666573579361080525984562);
    u = u * t + (0.0769219538311769618355029);
    u = u * t + (-0.090908995008245008229153);
    u = u * t + (0.111111105648261418443745);
    u = u * t + (-0.14285714266771329383765);
    u = u * t + (0.199999999996591265594148);
    u = u * t + (-0.333333333333311110369124);

    t = s + s * (t * u);

    if ((q & 1) != 0) t = 1.570796326794896557998982 - t;
    if ((q & 2) != 0) t = -t;

    return t;
}

__inline double2 atan2k_u1(double2 y, double2 x) {
    double u;
    double2 s, t;
    int q = 0;
    
    if (x.x < 0) { x.x = -x.x; x.y = -x.y; q = -2; }
    if (y.x > x.x) { t = x; x = y; y.x = -t.x; y.y = -t.y; q += 1; }
    
    s = dddiv_d2_d2_d2(y, x);
    t = ddsqu_d2_d2(s);
    t = ddnormalize_d2_d2(t);
    
    u = 1.06298484191448746607415e-05;
    u = mla(u, t.x, -0.000125620649967286867384336);
    u = mla(u, t.x, 0.00070557664296393412389774);
    u = mla(u, t.x, -0.00251865614498713360352999);
    u = mla(u, t.x, 0.00646262899036991172313504);
    u = mla(u, t.x, -0.0128281333663399031014274);
    u = mla(u, t.x, 0.0208024799924145797902497);
    u = mla(u, t.x, -0.0289002344784740315686289);
    u = mla(u, t.x, 0.0359785005035104590853656);
    u = mla(u, t.x, -0.041848579703592507506027);
    u = mla(u, t.x, 0.0470843011653283988193763);
    u = mla(u, t.x, -0.0524914210588448421068719);
    u = mla(u, t.x, 0.0587946590969581003860434);
    u = mla(u, t.x, -0.0666620884778795497194182);
    u = mla(u, t.x, 0.0769225330296203768654095);
    u = mla(u, t.x, -0.0909090442773387574781907);
    u = mla(u, t.x, 0.111111108376896236538123);
    u = mla(u, t.x, -0.142857142756268568062339);
    u = mla(u, t.x, 0.199999999997977351284817);
    u = mla(u, t.x, -0.333333333333317605173818);
    
    t = ddmul_d2_d2_d(t, u);
    t = ddmul_d2_d2_d2(s, ddadd_d2_d_d2(1, t));
    if (xfabs(s.x) < 1e-200) t = s;
    t = ddadd2_d2_d2_d2(ddmul_d2_d2_d(dd(1.570796326794896557998982, 6.12323399573676603586882e-17), q), t);
    
    return t;
}

__inline double xatan2_u1(double y, double x) {
    double2 d = atan2k_u1(dd(xfabs(y), 0), dd(x, 0));
    double r = d.x + d.y;
    
    r = mulsign(r, x);
    if (xisinf(x) || x == 0) r = rtengine::RT_PI_2 - (xisinf(x) ? (sign(x) * (rtengine::RT_PI_2)) : 0);
    if (xisinf(y)          ) r = rtengine::RT_PI_2 - (xisinf(x) ? (sign(x) * (rtengine::RT_PI*1/4)) : 0);
    if (             y == 0) r = (sign(x) == -1 ? rtengine::RT_PI : 0);
    
    return xisnan(x) || xisnan(y) ? rtengine::RT_NAN : mulsign(r, y);
}

__inline double xasin_u1(double d) {
    double2 d2 = atan2k_u1(dd(xfabs(d), 0), ddsqrt_d2_d2(ddmul_d2_d2_d2(ddadd_d2_d_d(1, d), ddadd_d2_d_d(1,-d))));
    double r = d2.x + d2.y;
    if (xfabs(d) == 1) r = 1.570796326794896557998982;
    return mulsign(r, d);
}

__inline double xacos_u1(double d) {
    double2 d2 = atan2k_u1(ddsqrt_d2_d2(ddmul_d2_d2_d2(ddadd_d2_d_d(1, d), ddadd_d2_d_d(1,-d))), dd(xfabs(d), 0));
    d2 = ddscale_d2_d2_d(d2, mulsign(1, d));
    if (xfabs(d) == 1) d2 = dd(0, 0);
    if (sign(d) == -1) d2 = ddadd_d2_d2_d2(dd(3.141592653589793116, 1.2246467991473532072e-16), d2);
    return d2.x + d2.y;
}

__inline double xatan_u1(double d) {
    double2 d2 = atan2k_u1(dd(xfabs(d), 0), dd(1, 0));
    double r = d2.x + d2.y;
    if (xisinf(d)) r = 1.570796326794896557998982;
    return mulsign(r, d);
}

__inline double xsin(double d) {
    double u, s, t = d;

    int qh = xtrunc(d * (rtengine::RT_1_PI / (1 << 24)));
    int ql = xrint(d * rtengine::RT_1_PI - qh * (double)(1 << 24));

    d = mla(qh, -PI_A * (1 << 24), d);
    d = mla(ql, -PI_A,             d);
    d = mla(qh, -PI_B * (1 << 24), d);
    d = mla(ql, -PI_B,             d);
    d = mla(qh, -PI_C * (1 << 24), d);
    d = mla(ql, -PI_C,             d);
    d = mla((double)qh * (1 << 24) + ql, -PI_D, d);

    s = d * d;

    if ((ql & 1) != 0) d = -d;

    u = -7.97255955009037868891952e-18;
    u = mla(u, s, 2.81009972710863200091251e-15);
    u = mla(u, s, -7.64712219118158833288484e-13);
    u = mla(u, s, 1.60590430605664501629054e-10);
    u = mla(u, s, -2.50521083763502045810755e-08);
    u = mla(u, s, 2.75573192239198747630416e-06);
    u = mla(u, s, -0.000198412698412696162806809);
    u = mla(u, s, 0.00833333333333332974823815);
    u = mla(u, s, -0.166666666666666657414808);

    u = mla(s, u * d, d);

    if (!xisinf(t) && (xisnegzero(t) || xfabs(t) > TRIGRANGEMAX)) u = -0.0;

    return u;
}

__inline double xsin_u1(double d) {
    double u;
    double2 s, t, x;

    int qh = xtrunc(d * (rtengine::RT_1_PI / (1 << 24)));
    int ql = xrint(d * rtengine::RT_1_PI - qh * (double)(1 << 24));

    s = ddadd2_d2_d_d (d, qh * (-PI_A * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_A            ));
    s = ddadd2_d2_d2_d(s, qh * (-PI_B * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_B            ));
    s = ddadd2_d2_d2_d(s, qh * (-PI_C * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_C            ));
    s = ddadd2_d2_d2_d(s, ((double)qh * (1 << 24) + ql) * -PI_D);

    t = s;
    s = ddsqu_d2_d2(s);

    u = 2.72052416138529567917983e-15;
    u = mla(u, s.x, -7.6429259411395447190023e-13);
    u = mla(u, s.x, 1.60589370117277896211623e-10);
    u = mla(u, s.x, -2.5052106814843123359368e-08);
    u = mla(u, s.x, 2.75573192104428224777379e-06);
    u = mla(u, s.x, -0.000198412698412046454654947);
    u = mla(u, s.x, 0.00833333333333318056201922);

    x = ddadd_d2_d_d2(1, ddmul_d2_d2_d2(ddadd_d2_d_d(-0.166666666666666657414808, u * s.x), s));

    x = ddmul_d2_d2_d2(t, x);
    u = x.x + x.y;

    if ((ql & 1) != 0) u = -u;
    if (!xisinf(d) && (xisnegzero(d) || xfabs(d) > TRIGRANGEMAX)) u = -0.0;

    return u;
}

__inline double xcos(double d) {
    double u, s, t = d;
    
    int qh = xtrunc(d * (rtengine::RT_1_PI / (1LL << 23)) - 0.5 * (rtengine::RT_1_PI / (1LL << 23)));
    int ql = 2*xrint(d * rtengine::RT_1_PI - 0.5 - qh * (double)(1LL << 23))+1;
    
    d = mla(qh, -PI_A*0.5*(1LL << 24), d);
    d = mla(ql, -PI_A*0.5,             d);
    d = mla(qh, -PI_B*0.5*(1LL << 24), d);
    d = mla(ql, -PI_B*0.5,             d);
    d = mla(qh, -PI_C*0.5*(1LL << 24), d);
    d = mla(ql, -PI_C*0.5,             d);
    d = mla((double)qh*(1LL << 24) + ql , -PI_D*0.5, d);
    
    s = d * d;
    
    if ((ql & 2) == 0) d = -d;
    
    u = -7.97255955009037868891952e-18;
    u = mla(u, s, 2.81009972710863200091251e-15);
    u = mla(u, s, -7.64712219118158833288484e-13);
    u = mla(u, s, 1.60590430605664501629054e-10);
    u = mla(u, s, -2.50521083763502045810755e-08);
    u = mla(u, s, 2.75573192239198747630416e-06);
    u = mla(u, s, -0.000198412698412696162806809);
    u = mla(u, s, 0.00833333333333332974823815);
    u = mla(u, s, -0.166666666666666657414808);
    
    u = mla(s, u * d, d);
    
    if (!xisinf(t) && xfabs(t) > TRIGRANGEMAX) u = 0.0;
    
    return u;
}

__inline double xcos_u1(double d) {
    double u;
    double2 s, t, x;
    
    d = xfabs(d);
    
    int qh = xtrunc(d * (rtengine::RT_1_PI / (1LL << (23))) - 0.5 * (rtengine::RT_1_PI / (1LL << (23))));
    int ql = 2*xrint(d * rtengine::RT_1_PI - 0.5 - qh * (double)(1LL << (23)))+1;
    
    s = ddadd2_d2_d_d (d, qh * (-PI_A*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_A*0.5            ));
    s = ddadd2_d2_d2_d(s, qh * (-PI_B*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_B*0.5            ));
    s = ddadd2_d2_d2_d(s, qh * (-PI_C*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_C*0.5            ));
    s = ddadd2_d2_d2_d(s, ((double)qh * (1 << 24) + ql) * (-PI_D*0.5));
    
    t = s;
    s = ddsqu_d2_d2(s);
    
    u = 2.72052416138529567917983e-15;
    u = mla(u, s.x, -7.6429259411395447190023e-13);
    u = mla(u, s.x, 1.60589370117277896211623e-10);
    u = mla(u, s.x, -2.5052106814843123359368e-08);
    u = mla(u, s.x, 2.75573192104428224777379e-06);
    u = mla(u, s.x, -0.000198412698412046454654947);
    u = mla(u, s.x, 0.00833333333333318056201922);
    
    x = ddadd_d2_d_d2(1, ddmul_d2_d2_d2(ddadd_d2_d_d(-0.166666666666666657414808, u * s.x), s));
    
    x = ddmul_d2_d2_d2(t, x);
    
    u = x.x + x.y;
    
    if ((((int)ql) & 2) == 0) u = -u;
    if (!xisinf(d) && d > TRIGRANGEMAX) u = 0.0;
    
    return u;
}

__inline double2 xsincos(double d) {
    double u, s, t;
    double2 r;
    
    s = d;
    
    int qh = xtrunc(d * ((2 * rtengine::RT_1_PI) / (1 << 24)));
    int ql = xrint(d * (2 * rtengine::RT_1_PI) - qh * (double)(1 << 24));
    
    s = mla(qh, -PI_A * 0.5 * (1 << 24), s);
    s = mla(ql, -PI_A * 0.5,             s);
    s = mla(qh, -PI_B * 0.5 * (1 << 24), s);
    s = mla(ql, -PI_B * 0.5,             s);
    s = mla(qh, -PI_C * 0.5 * (1 << 24), s);
    s = mla(ql, -PI_C * 0.5,             s);
    s = mla((double)qh * (1 << 24) + ql, -PI_D * 0.5, s);
    
    t = s;
    
    s = s * s;
    
    u = 1.58938307283228937328511e-10;
    u = mla(u, s, -2.50506943502539773349318e-08);
    u = mla(u, s, 2.75573131776846360512547e-06);
    u = mla(u, s, -0.000198412698278911770864914);
    u = mla(u, s, 0.0083333333333191845961746);
    u = mla(u, s, -0.166666666666666130709393);
    u = u * s * t;
    
    r.x = t + u;
    
    if (xisnegzero(d)) r.x = -0.0;
    
    u = -1.13615350239097429531523e-11;
    u = mla(u, s, 2.08757471207040055479366e-09);
    u = mla(u, s, -2.75573144028847567498567e-07);
    u = mla(u, s, 2.48015872890001867311915e-05);
    u = mla(u, s, -0.00138888888888714019282329);
    u = mla(u, s, 0.0416666666666665519592062);
    u = mla(u, s, -0.5);
    
    r.y = u * s + 1;
    
    if ((ql & 1) != 0) { s = r.y; r.y = r.x; r.x = s; }
    if ((ql & 2) != 0) { r.x = -r.x; }
    if (((ql+1) & 2) != 0) { r.y = -r.y; }
    
    if (xisinf(d)) { r.x = r.y = rtengine::RT_NAN; }
    if (!xisinf(d) && xfabs(d) > TRIGRANGEMAX) { r.x = r.y = 0; }
    
    return r;
}

__inline double2 xsincos_u1(double d) {
    double u;
    double2 r, s, t, x;
    
    int qh = xtrunc(d * ((2 * rtengine::RT_1_PI) / (1 << 24)));
    int ql = xrint(d * (2 * rtengine::RT_1_PI) - qh * (double)(1 << 24));
    
    s = ddadd2_d2_d_d (d, qh * (-PI_A*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_A*0.5            ));
    s = ddadd2_d2_d2_d(s, qh * (-PI_B*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_B*0.5            ));
    s = ddadd2_d2_d2_d(s, qh * (-PI_C*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_C*0.5            ));
    s = ddadd2_d2_d2_d(s, ((double)qh * (1 << 24) + ql) * (-PI_D*0.5));
    
    t = s;
    s = ddsqu_d2_d2(s);
    s.x = s.x + s.y;
    
    u = 1.58938307283228937328511e-10;
    u = mla(u, s.x, -2.50506943502539773349318e-08);
    u = mla(u, s.x, 2.75573131776846360512547e-06);
    u = mla(u, s.x, -0.000198412698278911770864914);
    u = mla(u, s.x, 0.0083333333333191845961746);
    u = mla(u, s.x, -0.166666666666666130709393);
    
    u *= s.x * t.x;
    
    x = ddadd_d2_d2_d(t, u);
    r.x = x.x + x.y;
    
    if (xisnegzero(d)) r.x = -0.0;
    
    u = -1.13615350239097429531523e-11;
    u = mla(u, s.x, 2.08757471207040055479366e-09);
    u = mla(u, s.x, -2.75573144028847567498567e-07);
    u = mla(u, s.x, 2.48015872890001867311915e-05);
    u = mla(u, s.x, -0.00138888888888714019282329);
    u = mla(u, s.x, 0.0416666666666665519592062);
    u = mla(u, s.x, -0.5);
    
    x = ddadd_d2_d_d2(1, ddmul_d2_d_d(s.x, u));
    r.y = x.x + x.y;
    
    if ((ql & 1) != 0) { u = r.y; r.y = r.x; r.x = u; }
    if ((ql & 2) != 0) { r.x = -r.x; }
    if (((ql+1) & 2) != 0) { r.y = -r.y; }
    
    if (xisinf(d)) { r.x = r.y = rtengine::RT_NAN; }
    if (!xisinf(d) && xfabs(d) > TRIGRANGEMAX) { r.x = r.y = 0; }
    
    return r;
}

__inline double xtan(double d) {
    double u, s, x;
    
    int qh = xtrunc(d * ((2 * rtengine::RT_1_PI) / (1 << 24)));
    int ql = xrint(d * (2 * rtengine::RT_1_PI) - qh * (double)(1 << 24));
    
    x = mla(qh, -PI_A * 0.5 * (1 << 24), d);
    x = mla(ql, -PI_A * 0.5,             x);
    x = mla(qh, -PI_B * 0.5 * (1 << 24), x);
    x = mla(ql, -PI_B * 0.5,             x);
    x = mla(qh, -PI_C * 0.5 * (1 << 24), x);
    x = mla(ql, -PI_C * 0.5,             x);
    x = mla((double)qh * (1 << 24) + ql, -PI_D * 0.5, x);
    
    s = x * x;
    
    if ((ql & 1) != 0) x = -x;
    
    u = 9.99583485362149960784268e-06;
    u = mla(u, s, -4.31184585467324750724175e-05);
    u = mla(u, s, 0.000103573238391744000389851);
    u = mla(u, s, -0.000137892809714281708733524);
    u = mla(u, s, 0.000157624358465342784274554);
    u = mla(u, s, -6.07500301486087879295969e-05);
    u = mla(u, s, 0.000148898734751616411290179);
    u = mla(u, s, 0.000219040550724571513561967);
    u = mla(u, s, 0.000595799595197098359744547);
    u = mla(u, s, 0.00145461240472358871965441);
    u = mla(u, s, 0.0035923150771440177410343);
    u = mla(u, s, 0.00886321546662684547901456);
    u = mla(u, s, 0.0218694899718446938985394);
    u = mla(u, s, 0.0539682539049961967903002);
    u = mla(u, s, 0.133333333334818976423364);
    u = mla(u, s, 0.333333333333320047664472);
    
    u = mla(s, u * x, x);
    
    if ((ql & 1) != 0) u = 1.0 / u;
    
    if (xisinf(d)) u = rtengine::RT_NAN;
    
    return u;
}

__inline double xtan_u1(double d) {
    double u;
    double2 s, t, x;
    
    int qh = xtrunc(d * (rtengine::RT_PI_2 / (1 << 24)));
    s = ddadd2_d2_d2_d(ddmul_d2_d2_d(dd(M_2_PI_H, M_2_PI_L), d), (d < 0 ? -0.5 : 0.5) - qh * (double)(1 << 24));
    int ql = s.x + s.y;
    
    s = ddadd2_d2_d_d (d, qh * (-PI_A*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_A*0.5            ));
    s = ddadd2_d2_d2_d(s, qh * (-PI_B*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_B*0.5            ));
    s = ddadd2_d2_d2_d(s, qh * (-PI_C*0.5 * (1 << 24)));
    s = ddadd2_d2_d2_d(s, ql * (-PI_C*0.5            ));
    s = ddadd2_d2_d2_d(s, ((double)qh * (1 << 24) + ql) * (-PI_D*0.5));
    
    if ((ql & 1) != 0) s = ddneg_d2_d2(s);
    
    t = s;
    s = ddsqu_d2_d2(s);
    
    u = 1.01419718511083373224408e-05;
    u = mla(u, s.x, -2.59519791585924697698614e-05);
    u = mla(u, s.x, 5.23388081915899855325186e-05);
    u = mla(u, s.x, -3.05033014433946488225616e-05);
    u = mla(u, s.x, 7.14707504084242744267497e-05);
    u = mla(u, s.x, 8.09674518280159187045078e-05);
    u = mla(u, s.x, 0.000244884931879331847054404);
    u = mla(u, s.x, 0.000588505168743587154904506);
    u = mla(u, s.x, 0.00145612788922812427978848);
    u = mla(u, s.x, 0.00359208743836906619142924);
    u = mla(u, s.x, 0.00886323944362401618113356);
    u = mla(u, s.x, 0.0218694882853846389592078);
    u = mla(u, s.x, 0.0539682539781298417636002);
    u = mla(u, s.x, 0.133333333333125941821962);
    
    x = ddadd_d2_d_d2(1, ddmul_d2_d2_d2(ddadd_d2_d_d(0.333333333333334980164153, u * s.x), s));
    x = ddmul_d2_d2_d2(t, x);
    
    if ((ql & 1) != 0) x = ddrec_d2_d2(x);
    
    u = x.x + x.y;
    
    if (!xisinf(d) && (xisnegzero(d) || xfabs(d) > TRIGRANGEMAX)) u = -0.0;
    
    return u;
}

__inline double xlog(double d) {
    double x, x2, t, m;
    int e;

    e = ilogbk(d * (1.0/0.75));
    m = ldexpk(d, -e);

    x = (m-1) / (m+1);
    x2 = x * x;

    t = 0.153487338491425068243146;
    t = mla(t, x2, 0.152519917006351951593857);
    t = mla(t, x2, 0.181863266251982985677316);
    t = mla(t, x2, 0.222221366518767365905163);
    t = mla(t, x2, 0.285714294746548025383248);
    t = mla(t, x2, 0.399999999950799600689777);
    t = mla(t, x2, 0.6666666666667778740063);
    t = mla(t, x2, 2);
    
    x = x * t + 0.693147180559945286226764 * e;

    if (xispinf(d)) x = rtengine::RT_INFINITY;
    if (d < 0) x = rtengine::RT_NAN;
    if (d == 0) x = -rtengine::RT_INFINITY;

    return x;
}

__inline double xexp(double d) {
    int q = (int)xrint(d * R_LN2);
    double s, u;

    s = mla(q, -L2U, d);
    s = mla(q, -L2L, s);

    u = 2.08860621107283687536341e-09;
    u = mla(u, s, 2.51112930892876518610661e-08);
    u = mla(u, s, 2.75573911234900471893338e-07);
    u = mla(u, s, 2.75572362911928827629423e-06);
    u = mla(u, s, 2.4801587159235472998791e-05);
    u = mla(u, s, 0.000198412698960509205564975);
    u = mla(u, s, 0.00138888888889774492207962);
    u = mla(u, s, 0.00833333333331652721664984);
    u = mla(u, s, 0.0416666666666665047591422);
    u = mla(u, s, 0.166666666666666851703837);
    u = mla(u, s, 0.5);

    u = s * s * u + s + 1;
    u = ldexpk(u, q);

    if (d < -1000) u = 0;

    return u;
}

__inline double2 logk(double d) {
    double2 x, x2;
    double m, t;
    int e;

    e = ilogbk(d * (1.0/0.75));
    m = ldexpk(d, -e);

    x = dddiv_d2_d2_d2(ddadd2_d2_d_d(-1, m), ddadd2_d2_d_d(1, m));
    x2 = ddsqu_d2_d2(x);

    t = 0.116255524079935043668677;
    t = mla(t, x2.x, 0.103239680901072952701192);
    t = mla(t, x2.x, 0.117754809412463995466069);
    t = mla(t, x2.x, 0.13332981086846273921509);
    t = mla(t, x2.x, 0.153846227114512262845736);
    t = mla(t, x2.x, 0.181818180850050775676507);
    t = mla(t, x2.x, 0.222222222230083560345903);
    t = mla(t, x2.x, 0.285714285714249172087875);
    t = mla(t, x2.x, 0.400000000000000077715612);
    double2 c = dd(0.666666666666666629659233, 3.80554962542412056336616e-17);
    
    return ddadd2_d2_d2_d2(ddmul_d2_d2_d(dd(0.693147180559945286226764, 2.319046813846299558417771e-17), e),
            ddadd2_d2_d2_d2(ddscale_d2_d2_d(x, 2),
                    ddmul_d2_d2_d2(ddmul_d2_d2_d2(x2, x),
                            ddadd2_d2_d2_d2(ddmul_d2_d2_d(x2, t), c))));
}

__inline double xlog_u1(double d) {
    double2 s = logk(d);
    double x = s.x + s.y;
    
    if (xisinf(d)) x = rtengine::RT_INFINITY;
    if (d < 0) x = rtengine::RT_NAN;
    if (d == 0) x = -rtengine::RT_INFINITY;
    
    return x;
}

__inline double expk(double2 d) {
    int q = (int)rint((d.x + d.y) * R_LN2);
    double2 s, t;
    double u;

    s = ddadd2_d2_d2_d(d, q * -L2U);
    s = ddadd2_d2_d2_d(s, q * -L2L);
    
    s = ddnormalize_d2_d2(s);
    
    u = 2.51069683420950419527139e-08;
    u = mla(u, s.x, 2.76286166770270649116855e-07);
    u = mla(u, s.x, 2.75572496725023574143864e-06);
    u = mla(u, s.x, 2.48014973989819794114153e-05);
    u = mla(u, s.x, 0.000198412698809069797676111);
    u = mla(u, s.x, 0.0013888888939977128960529);
    u = mla(u, s.x, 0.00833333333332371417601081);
    u = mla(u, s.x, 0.0416666666665409524128449);
    u = mla(u, s.x, 0.166666666666666740681535);
    u = mla(u, s.x, 0.500000000000000999200722);
    
    t = ddadd_d2_d2_d2(s, ddmul_d2_d2_d(ddsqu_d2_d2(s), u));
    
    t = ddadd_d2_d_d2(1, t);
    
    u = ldexpk(t.x + t.y, q);
    
    if (d.x < -1000) u = 0;
    
    return u;
}

__inline double xpow(double x, double y) {
    int yisint = xisint(y);
    int yisodd = (1 & (int)y) != 0 && yisint;

    double result = expk(ddmul_d2_d2_d(logk(xfabs(x)), y));

    result = xisnan(result) ? rtengine::RT_INFINITY : result;
    result *= (x > 0 ? 1 : (!yisint ? rtengine::RT_NAN : (yisodd ? -1 : 1)));

    double efx = mulsign(xfabs(x) - 1, y);
    if (xisinf(y)) result = efx < 0 ? 0.0 : (efx == 0 ? 1.0 : rtengine::RT_INFINITY);
    if (xisinf(x) || x == 0) result = (yisodd ? sign(x) : 1) * ((x == 0 ? -y : y) < 0 ? 0 : rtengine::RT_INFINITY);
    if (xisnan(x) || xisnan(y)) result = rtengine::RT_NAN;
    if (y == 0 || x == 1) result = 1;

    return result;
}

__inline double2 expk2(double2 d) {
    int q = (int)xrint((d.x + d.y) * R_LN2);
    double2 s, t;
    double u;

    s = ddadd2_d2_d2_d(d, q * -L2U);
    s = ddadd2_d2_d2_d(s, q * -L2L);
    
    u = 2.51069683420950419527139e-08;
    u = mla(u, s.x, 2.76286166770270649116855e-07);
    u = mla(u, s.x, 2.75572496725023574143864e-06);
    u = mla(u, s.x, 2.48014973989819794114153e-05);
    u = mla(u, s.x, 0.000198412698809069797676111);
    u = mla(u, s.x, 0.0013888888939977128960529);
    u = mla(u, s.x, 0.00833333333332371417601081);
    u = mla(u, s.x, 0.0416666666665409524128449);
    u = mla(u, s.x, 0.166666666666666740681535);
    u = mla(u, s.x, 0.500000000000000999200722);
    
    t = ddadd_d2_d2_d2(s, ddmul_d2_d2_d(ddsqu_d2_d2(s), u));
    
    t = ddadd_d2_d_d2(1, t);
    return ddscale_d2_d2_d(ddscale_d2_d2_d(t, 2), pow2i(q-1));
}

__inline double xsinh(double x) {
    double y = xfabs(x);
    double2 d = expk2(dd(y, 0));
    d = ddsub_d2_d2_d2(d, ddrec_d2_d2(d));
    y = (d.x + d.y) * 0.5;

    y = xfabs(x) > 710 ? rtengine::RT_INFINITY : y;
    y = xisnan(y) ? rtengine::RT_INFINITY : y;
    y = mulsign(y, x);
    y = xisnan(x) ? rtengine::RT_NAN : y;

    return y;
}

__inline double xcosh(double x) {
    double y = xfabs(x);
    double2 d = expk2(dd(x, 0));
    d = ddadd_d2_d2_d2(d, ddrec_d2_d2(d));
    y = (d.x + d.y) * 0.5;

    y = xfabs(x) > 710 ? rtengine::RT_INFINITY : y;
    y = xisnan(y) ? rtengine::RT_INFINITY : y;
    y = xisnan(x) ? rtengine::RT_NAN : y;

    return y;
}

__inline double xtanh(double x) {
    double y = xfabs(x);
    double2 d = expk2(dd(y, 0));
    double2 e = ddrec_d2_d2(d);
    d = dddiv_d2_d2_d2(ddsub_d2_d2_d2(d, e), ddadd_d2_d2_d2(d, e));
    y = d.x + d.y;

    y = xfabs(x) > 18.714973875 ? 1.0 : y;
    y = xisnan(y) ? 1.0 : y;
    y = mulsign(y, x);
    y = xisnan(x) ? rtengine::RT_NAN : y;

    return y;
}

__inline double2 logk2(double2 d) {
    double2 x, x2, m;
    double t;
    int e;

    e = ilogbk(d.x * (1.0/0.75));
    m = ddscale_d2_d2_d(d, pow2i(-e));

    x = dddiv_d2_d2_d2(ddadd2_d2_d2_d(m, -1), ddadd2_d2_d2_d(m, 1));
    x2 = ddsqu_d2_d2(x);
    
    t = 0.13860436390467167910856;
    t = mla(t, x2.x, 0.131699838841615374240845);
    t = mla(t, x2.x, 0.153914168346271945653214);
    t = mla(t, x2.x, 0.181816523941564611721589);
    t = mla(t, x2.x, 0.22222224632662035403996);
    t = mla(t, x2.x, 0.285714285511134091777308);
    t = mla(t, x2.x, 0.400000000000914013309483);
    t = mla(t, x2.x, 0.666666666666664853302393);
    
    return ddadd2_d2_d2_d2(ddmul_d2_d2_d(dd(0.693147180559945286226764, 2.319046813846299558417771e-17), e),
            ddadd2_d2_d2_d2(ddscale_d2_d2_d(x, 2), ddmul_d2_d2_d(ddmul_d2_d2_d2(x2, x), t)));
}

__inline double xasinh(double x) {
    double y = xfabs(x);
    double2 d;
    
    d = y > 1 ? ddrec_d2_d(x) : dd(y, 0);
    d = ddsqrt_d2_d2(ddadd2_d2_d2_d(ddsqu_d2_d2(d), 1));
    d = y > 1 ? ddmul_d2_d2_d(d, y) : d;
    
    d = logk2(ddnormalize_d2_d2(ddadd_d2_d2_d(d, x)));
    y = d.x + d.y;
    
    y = (xfabs(x) > SQRT_DBL_MAX || xisnan(y)) ? mulsign(rtengine::RT_INFINITY, x) : y;
    y = xisnan(x) ? rtengine::RT_NAN : y;
    y = xisnegzero(x) ? -0.0 : y;
    
    return y;
}

__inline double xacosh(double x) {
    double2 d = logk2(ddadd2_d2_d2_d(ddmul_d2_d2_d2(ddsqrt_d2_d2(ddadd2_d2_d_d(x, 1)), ddsqrt_d2_d2(ddadd2_d2_d_d(x, -1))), x));
    double y = d.x + d.y;
    
    y = (x > SQRT_DBL_MAX || xisnan(y)) ? rtengine::RT_INFINITY : y;
    y = x == 1.0 ? 0.0 : y;
    y = x < 1.0 ? rtengine::RT_NAN : y;
    y = xisnan(x) ? rtengine::RT_NAN : y;
    
    return y;
}

__inline double xatanh(double x) {
    double y = xfabs(x);
    double2 d = logk2(dddiv_d2_d2_d2(ddadd2_d2_d_d(1, y), ddadd2_d2_d_d(1, -y)));
    y = y > 1.0 ? rtengine::RT_NAN : (y == 1.0 ? rtengine::RT_INFINITY : (d.x + d.y) * 0.5);

    y = xisinf(x) || xisnan(y) ? rtengine::RT_NAN : y;
    y = mulsign(y, x);
    y = xisnan(x) ? rtengine::RT_NAN : y;

    return y;
}

//

__inline double xfma(double x, double y, double z) {
    union {
        double f;
        long long int i;
    } tmp;

    tmp.f = x;
    tmp.i = (tmp.i + 0x4000000) & 0xfffffffff8000000LL;
    double xh = tmp.f, xl = x - xh;

    tmp.f = y;
    tmp.i = (tmp.i + 0x4000000) & 0xfffffffff8000000LL;
    double yh = tmp.f, yl = y - yh;

    double h = x * y;
    double l = xh * yh - h + xl * yh + xh * yl + xl * yl;

    double h2, l2, v;

    h2 = h + z;
    v = h2 - h;
    l2 = (h - (h2 - v)) + (z - v) + l;

    return h2 + l2;
}

__inline double xsqrt(double d) { // max error : 0.5 ulp
    double q = 1;

    if (d < 8.636168555094445E-78) {
        d *= 1.157920892373162E77;
        q = 2.9387358770557188E-39;
    }

    // http://en.wikipedia.org/wiki/Fast_inverse_square_root
    double x = longBitsToDouble(0x5fe6ec85e7de30da - (doubleToRawLongBits(d + 1e-320) >> 1));

    x = x * (1.5 - 0.5 * d * x * x);
    x = x * (1.5 - 0.5 * d * x * x);
    x = x * (1.5 - 0.5 * d * x * x);

    // You can change xfma to fma if fma is correctly implemented
    x = xfma(d * x, d * x, -d) * (x * -0.5) + d * x;

    return d == rtengine::RT_INFINITY ? rtengine::RT_INFINITY : x * q;
}

__inline double xcbrt(double d) { // max error : 2 ulps
    double x, y, q = 1.0;
    int e, r;

    e = ilogbk(xfabs(d))+1;
    d = ldexpk(d, -e);
    r = (e + 6144) % 3;
    q = (r == 1) ? 1.2599210498948731647672106 : q;
    q = (r == 2) ? 1.5874010519681994747517056 : q;
    q = ldexpk(q, (e + 6144) / 3 - 2048);

    q = mulsign(q, d);
    d = xfabs(d);

    x = -0.640245898480692909870982;
    x = x * d + 2.96155103020039511818595;
    x = x * d + -5.73353060922947843636166;
    x = x * d + 6.03990368989458747961407;
    x = x * d + -3.85841935510444988821632;
    x = x * d + 2.2307275302496609725722;

    y = x * x; y = y * y; x -= (d * y - x) * (1.0 / 3.0);
    y = d * x * x;
    y = (y - (2.0 / 3.0) * y * (y * x - 1)) * q;

    return y;
}

__inline double xcbrt_u1(double d) {
    double x, y, z;
    double2 q2 = dd(1, 0), u, v;
    int e, r;
    
    e = ilogbk(xfabs(d))+1;
    d = ldexpk(d, -e);
    r = (e + 6144) % 3;
    q2 = (r == 1) ? dd(1.2599210498948731907, -2.5899333753005069177e-17) : q2;
    q2 = (r == 2) ? dd(1.5874010519681995834, -1.0869008194197822986e-16) : q2;
    
    q2.x = mulsign(q2.x, d); q2.y = mulsign(q2.y, d);
    d = xfabs(d);
    
    x = -0.640245898480692909870982;
    x = x * d + 2.96155103020039511818595;
    x = x * d + -5.73353060922947843636166;
    x = x * d + 6.03990368989458747961407;
    x = x * d + -3.85841935510444988821632;
    x = x * d + 2.2307275302496609725722;
    
    y = x * x; y = y * y; x -= (d * y - x) * (1.0 / 3.0);
    
    z = x;
    
    u = ddmul_d2_d_d(x, x);
    u = ddmul_d2_d2_d2(u, u);
    u = ddmul_d2_d2_d(u, d);
    u = ddadd2_d2_d2_d(u, -x);
    y = u.x + u.y;
    
    y = -2.0 / 3.0 * y * z;
    v = ddadd2_d2_d2_d(ddmul_d2_d_d(z, z), y);
    v = ddmul_d2_d2_d(v, d);
    v = ddmul_d2_d2_d2(v, q2);
    z = ldexp(v.x + v.y, (e + 6144) / 3 - 2048);
    
    if (xisinf(d)) { z = mulsign(rtengine::RT_INFINITY, q2.x); }
    if (d == 0) { z = mulsign(0, q2.x); }
    
    return z;
}

__inline double xexp2(double a) {
    double u = expk(ddmul_d2_d2_d(dd(0.69314718055994528623, 2.3190468138462995584e-17), a));
    if (a > 1024) u = rtengine::RT_INFINITY; // log2(DBL_MAX)
    if (xisminf(a)) u = 0;
    return u;
}

__inline double xexp10(double a) {
    double u = expk(ddmul_d2_d2_d(dd(2.3025850929940459011, -2.1707562233822493508e-16), a));
    if (a > 308.254715559916743850652254) u = rtengine::RT_INFINITY; // log10(DBL_MAX)
    if (xisminf(a)) u = 0;
    return u;
}

__inline double xexpm1(double a) {
    double2 d = ddadd2_d2_d2_d(expk2(dd(a, 0)), -1.0);
    double x = d.x + d.y;
    if (a > 709.782712893383996732223) x = rtengine::RT_INFINITY; // log(DBL_MAX)
    if (a < -36.736800569677101399113302437) x = -1; // log(1 - nexttoward(1, 0))
    if (xisnegzero(a)) x = -0.0;
    return x;
}

__inline double xlog10(double a) {
    double2 d = ddmul_d2_d2_d2(logk(a), dd(0.43429448190325176116, 6.6494347733425473126e-17));
    double x = d.x + d.y;

    if (xisinf(a)) x = rtengine::RT_INFINITY;
    if (a < 0) x = rtengine::RT_NAN;
    if (a == 0) x = -rtengine::RT_INFINITY;

    return x;
}

__inline double xlog1p(double a) {
    double2 d = logk2(ddadd2_d2_d_d(a, 1));
    double x = d.x + d.y;

    if (a > 1e+307) x = rtengine::RT_INFINITY;
    if (a < -1) x = rtengine::RT_NAN;
    if (a == -1) x = -rtengine::RT_INFINITY;
    if (xisnegzero(a)) x = -0.0;

    return x;
}

//
// Single precision functions (see sleef-2.121/purec/sleefsp.c)
// Some modifications made to embed SSE versions
// Not all functions included, only those needed
//

#define PI_Af 3.140625f
#define PI_Bf 0.0009670257568359375f
#define PI_Cf 6.2771141529083251953e-07f
#define PI_Df 1.2154201256553420762e-10f

#define TRIGRANGEMAXf 1e+7 // 39000

#define L2Uf 0.693145751953125f
#define L2Lf 1.428606765330187045e-06f

#define R_LN2f 1.442695040888963407359924681001892137426645954152985934135449406931f

__inline int32_t floatToRawIntBits(float d) {
    union {
        float f;
        int32_t i;
    } tmp;
    tmp.f = d;
    return tmp.i;
}

__inline float intBitsToFloat(int32_t i) {
    union {
        float f;
        int32_t i;
    } tmp;
    tmp.i = i;
    return tmp.f;
}

__inline float xfabsf(float x) {
    return intBitsToFloat(0x7fffffffL & floatToRawIntBits(x));
}

__inline float mulsignf(float x, float y) {
    return intBitsToFloat(floatToRawIntBits(x) ^ (floatToRawIntBits(y) & (1 << 31)));
}

__inline float signf(float d) { return mulsignf(1, d); }
__inline float mlaf(float x, float y, float z) { return x * y + z; }

// RT modification
#ifdef __SSE2__
__inline int xrintf(float x) {
    return _mm_cvt_ss2si(_mm_set_ss(x));
}
#else
__inline int xrintf(float x) {
    return x + (x < 0 ? -0.5f : 0.5f);
}
#endif

__inline int xisnanf(float x) { return x != x; }
__inline int xisinff(float x) { return x == rtengine::RT_INFINITY_F || x == -rtengine::RT_INFINITY_F; }
__inline int xisminff(float x) { return x == -rtengine::RT_INFINITY_F; }
__inline int xispinff(float x) { return x == rtengine::RT_INFINITY_F; }
__inline int xisnegzerof(float x) { return floatToRawIntBits(x) == floatToRawIntBits(-0.0); }

__inline int ilogbkf(float d) {
    int m = d < 5.421010862427522E-20f;
    d = m ? 1.8446744073709552E19f * d : d;
    int q = (floatToRawIntBits(d) >> 23) & 0xff;
    q = m ? q - (64 + 0x7f) : q - 0x7f;
    return q;
}

__inline float ldexpkf(float x, int q) {
    float u;
    int m;
    m = q >> 31;
    m = (((m + q) >> 6) - m) << 4;
    q = q - (m << 2);
    m += 127;
    m = m <   0 ?   0 : m;
    m = m > 255 ? 255 : m;
    u = intBitsToFloat(((int32_t)m) << 23);
    x = x * u * u * u * u;
    u = intBitsToFloat(((int32_t)(q + 0x7f)) << 23);
    return x * u;
}

__inline float xcbrtf(float d) {
    float x, y, q = 1.0f;
    int e, r;

    e = ilogbkf(xfabsf(d))+1;
    d = ldexpkf(d, -e);
    r = (e + 6144) % 3;
    q = (r == 1) ? 1.2599210498948731647672106f : q;
    q = (r == 2) ? 1.5874010519681994747517056f : q;
    q = ldexpkf(q, (e + 6144) / 3 - 2048);

    q = mulsignf(q, d);
    d = xfabsf(d);

    x = -0.601564466953277587890625f;
    x = mlaf(x, d, 2.8208892345428466796875f);
    x = mlaf(x, d, -5.532182216644287109375f);
    x = mlaf(x, d, 5.898262500762939453125f);
    x = mlaf(x, d, -3.8095417022705078125f);
    x = mlaf(x, d, 2.2241256237030029296875f);

    y = d * x * x;
    y = (y - (2.0f / 3.0f) * y * (y * x - 1.0f)) * q;

    return y;
}

__inline float xsinf(float d) {
    int q;
    float u, s, t = d;

    q = xrintf(d * rtengine::RT_1_PI_F);

    d = mlaf(q, -PI_Af, d);
    d = mlaf(q, -PI_Bf, d);
    d = mlaf(q, -PI_Cf, d);
    d = mlaf(q, -PI_Df, d);

    if (floatToRawIntBits(d) == floatToRawIntBits(-0.0f)) s = -0.0f;
    if ((q & 1) != 0) d = -d;

    u = 2.6083159809786593541503e-06f;
    u = mlaf(u, s, -0.0001981069071916863322258f);
    u = mlaf(u, s, 0.00833307858556509017944336f);
    u = mlaf(u, s, -0.166666597127914428710938f);

    u = mlaf(s, u * d, d);
    
    if (xisnegzerof(t) || xfabsf(t) > TRIGRANGEMAXf) u = -0.0f;
    if (xisinff(t)) u = rtengine::RT_NAN_F;

    return u;
}

__inline float xcosf(float d) {
#ifdef __SSE2__
    // faster than scalar version
    return xcosf(_mm_set_ss(d))[0];
#else
    int q;
    float u, s ,t = d;

    q = 1 + 2*xrintf(d * rtengine::RT_1_PI_F - 0.5f);

    d = mlaf(q, -PI_Af*0.5f, d);
    d = mlaf(q, -PI_Bf*0.5f, d);
    d = mlaf(q, -PI_Cf*0.5f, d);
    d = mlaf(q, -PI_Df*0.5f, d);

    s = d * d;

    if ((q & 2) == 0) d = -d;

    u = 2.6083159809786593541503e-06f;
    u = mlaf(u, s, -0.0001981069071916863322258f);
    u = mlaf(u, s, 0.00833307858556509017944336f);
    u = mlaf(u, s, -0.166666597127914428710938f);

    u = mlaf(s, u * d, d);
    
    if (xfabsf(t) > TRIGRANGEMAXf) u = 0.0f;
    if (xisinff(t)) u = rtengine::RT_NAN_F;

    return u;
#endif
}

__inline float2 xsincosf(float d) {
#ifdef __SSE2__
    // faster than scalar version
    vfloat2 res = xsincosf(_mm_set_ss(d));
    return {res.x[0], res.y[0]};
#else
    int q;
    float u, s, t;
    float2 r;

    q = xrintf(d * rtengine::RT_2_PI_F);

    s = d;

    s = mlaf(q, -PI_Af*0.5f, s);
    s = mlaf(q, -PI_Bf*0.5f, s);
    s = mlaf(q, -PI_Cf*0.5f, s);
    s = mlaf(q, -PI_Df*0.5f, s);

    t = s;

    s = s * s;

    u = -0.000195169282960705459117889f;
    u = mlaf(u, s, 0.00833215750753879547119141f);
    u = mlaf(u, s, -0.166666537523269653320312f);
    u = u * s * t;

    r.x = t + u;
    
    if (xisnegzerof(d)) r.x = -0.0f;

    u = -2.71811842367242206819355e-07f;
    u = mlaf(u, s, 2.47990446951007470488548e-05f);
    u = mlaf(u, s, -0.00138888787478208541870117f);
    u = mlaf(u, s, 0.0416666641831398010253906f);
    u = mlaf(u, s, -0.5f);

    r.y = u * s + 1;

    if ((q & 1) != 0) { s = r.y; r.y = r.x; r.x = s; }
    if ((q & 2) != 0) { r.x = -r.x; }
    if (((q+1) & 2) != 0) { r.y = -r.y; }

    if (xfabsf(d) > TRIGRANGEMAXf) { r.x = r.y = 0; }
    if (xisinff(d)) { r.x = r.y = rtengine::RT_NAN_F; }

    return r;
#endif
}

__inline float atan2kf(float y, float x) {
    float s, t, u;
    int q = 0;

    if (x < 0) { x = -x; q = -2; }
    if (y > x) { t = x; x = y; y = -t; q += 1; }

    s = y / x;
    t = s * s;

    u = 0.00282363896258175373077393f;
    u = mlaf(u, t, -0.0159569028764963150024414f);
    u = mlaf(u, t, 0.0425049886107444763183594f);
    u = mlaf(u, t, -0.0748900920152664184570312f);
    u = mlaf(u, t, 0.106347933411598205566406f);
    u = mlaf(u, t, -0.142027363181114196777344f);
    u = mlaf(u, t, 0.199926957488059997558594f);
    u = mlaf(u, t, -0.333331018686294555664062f);

    t = u * t * s + s;
    t = q * rtengine::RT_PI_F_2 + t;

    return t;
}

__inline float xatan2f(float y, float x) {
    float r = atan2kf(xfabsf(y), x);

    r = mulsignf(r, x);
    if (xisinff(x) || x == 0) r = rtengine::RT_PI_F/2 - (xisinff(x) ? (signf(x) * (float)(rtengine::RT_PI_F*.5f)) : 0);
    if (xisinff(y)          ) r = rtengine::RT_PI_F/2 - (xisinff(x) ? (signf(x) * (float)(rtengine::RT_PI_F*.25f)) : 0);
    if (              y == 0) r = (signf(x) == -1 ? rtengine::RT_PI_F : 0);

    return xisnanf(x) || xisnanf(y) ? rtengine::RT_NAN_F : mulsignf(r, y);
}

__inline float xlogf(float d) {
    float x, x2, t, m;
    int e;

    e = ilogbkf(d * (1.0f/0.75f));
    m = ldexpkf(d, -e);

    x = (m-1.0f) / (m+1.0f);
    x2 = x * x;

    t = 0.2392828464508056640625f;
    t = mlaf(t, x2, 0.28518211841583251953125f);
    t = mlaf(t, x2, 0.400005877017974853515625f);
    t = mlaf(t, x2, 0.666666686534881591796875f);
    t = mlaf(t, x2, 2.0f);

    x = x * t + 0.693147180559945286226764f * e;

    if (xispinff(d)) x = rtengine::RT_INFINITY_F;
    if (d < 0) x = rtengine::RT_NAN_F;
    if (d == 0) x = -rtengine::RT_INFINITY_F;

    return x;
}

// RT custom function used only in iptransform.cc
__inline float xlogf1(float d) { // does xlogf(vmaxf(d, 1.f)) but faster
    float x, x2, t, m;
    int e;

    e = ilogbkf(d * (1.0f/0.75f));
    m = ldexpkf(d, -e);

    x = (m-1.0f) / (m+1.0f);
    x2 = x * x;

    t = 0.2392828464508056640625f;
    t = mlaf(t, x2, 0.28518211841583251953125f);
    t = mlaf(t, x2, 0.400005877017974853515625f);
    t = mlaf(t, x2, 0.666666686534881591796875f);
    t = mlaf(t, x2, 2.0f);

    x = x * t + 0.693147180559945286226764f * e;

    if (xispinff(d)) x = rtengine::RT_INFINITY_F;
    if (d <= 1.f) x = 0;

    return x;
}

__inline float xexpf(float d) {
    int q = (int)xrintf(d * R_LN2f);
    float s, u;

    s = mlaf(q, -L2Uf, d);
    s = mlaf(q, -L2Lf, s);

    u = 0.000198527617612853646278381;
    u = mlaf(u, s, 0.00139304355252534151077271);
    u = mlaf(u, s, 0.00833336077630519866943359);
    u = mlaf(u, s, 0.0416664853692054748535156);
    u = mlaf(u, s, 0.166666671633720397949219);
    u = mlaf(u, s, 0.5);

    u = s * s * u + s + 1.0f;
    u = ldexpkf(u, q);

    if (d < -104) u = 0;

    return u;
}

//
// Functions not present in sleef, added by RT
//

__inline float xmul2f(float d) {
    union {
        float floatval;
        int intval;
    } uflint;
    uflint.floatval = d;
    if (uflint.intval & 0x7FFFFFFF) { // if f==0 do nothing
        uflint.intval += 1 << 23; // add 1 to the exponent
    }
    return uflint.floatval;
}

__inline float xdiv2f(float d) {
    union {
        float floatval;
        int intval;
    } uflint;
    uflint.floatval = d;
    if (uflint.intval & 0x7FFFFFFF) { // if f==0 do nothing
        uflint.intval -= 1 << 23; // sub 1 from the exponent
    }
    return uflint.floatval;
}

__inline float xdivf( float d, int n){
    union {
        float floatval;
        int intval;
    } uflint;
    uflint.floatval = d;
    if (uflint.intval & 0x7FFFFFFF) { // if f==0 do nothing
        uflint.intval -= n << 23; // add n to the exponent
    }
    return uflint.floatval;
}

__inline float xlin2log(float x, float base)
{
    constexpr float one(1);
    return xlogf(x * (base - one) + one) / xlogf(base);
}

__inline float xlog2lin(float x, float base)
{
    constexpr float one(1);
    return (xexpf(base * xlogf(x)) - one) / (base - one); 
}
