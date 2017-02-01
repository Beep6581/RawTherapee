////////////////////////////////////////////////////////////////
//
//  this code was taken from http://shibatch.sourceforge.net/
//  Many thanks to the author of original version: Naoki Shibata
//
//  This version contains modifications made by Ingo Weyrich
//
////////////////////////////////////////////////////////////////


#ifndef SLEEFSSEAVX
#define SLEEFSSEAVX

#include <assert.h>
//#include <math.h>
//#include <bits/nan.h>
//#include <bits/inf.h>
//#include "sleefsseavx.h"
#include "rt_math.h"
#ifdef __SSE2__
#include "helpersse2.h"

#ifdef ENABLE_AVX
#include "helperavx.h"
#endif

#ifdef __GNUC__
#define INLINE __inline
#else
#define INLINE inline
#endif

//

#define PI4_A .7853981554508209228515625
#define PI4_B .794662735614792836713604629039764404296875e-8
#define PI4_C .306161699786838294306516483068750264552437361480769e-16
#define M_4_PI 1.273239544735162542821171882678754627704620361328125

#define L2U .69314718055966295651160180568695068359375
#define L2L .28235290563031577122588448175013436025525412068e-12
#define R_LN2 1.442695040888963407359924681001892137426645954152985934135449406931

//

#define PI4_Af 0.78515625f
#define PI4_Bf 0.00024127960205078125f
#define PI4_Cf 6.3329935073852539062e-07f
#define PI4_Df 4.9604681473525147339e-10f

#define L2Uf 0.693145751953125f
#define L2Lf 1.428606765330187045e-06f
#define R_LN2f 1.442695040888963407359924681001892137426645954152985934135449406931f

#define INFINITYf ((float)rtengine::RT_INFINITY)
#define NANf ((float)rtengine::RT_NAN)

//

static INLINE vdouble vadd3(vdouble v0, vdouble v1, vdouble v2) {
  return vadd(vadd(v0, v1), v2);
}

static INLINE vdouble vadd4(vdouble v0, vdouble v1, vdouble v2, vdouble v3) {
  return vadd3(vadd(v0, v1), v2, v3);
}

static INLINE vdouble vadd5(vdouble v0, vdouble v1, vdouble v2, vdouble v3, vdouble v4) {
  return vadd4(vadd(v0, v1), v2, v3, v4);
}

static INLINE vdouble vadd6(vdouble v0, vdouble v1, vdouble v2, vdouble v3, vdouble v4, vdouble v5) {
  return vadd5(vadd(v0, v1), v2, v3, v4, v5);
}

static INLINE vdouble vadd7(vdouble v0, vdouble v1, vdouble v2, vdouble v3, vdouble v4, vdouble v5, vdouble v6) {
  return vadd6(vadd(v0, v1), v2, v3, v4, v5, v6);
}

static INLINE vdouble vsub3(vdouble v0, vdouble v1, vdouble v2) {
  return vsub(vsub(v0, v1), v2);
}

static INLINE vdouble vsub4(vdouble v0, vdouble v1, vdouble v2, vdouble v3) {
  return vsub3(vsub(v0, v1), v2, v3);
}

static INLINE vdouble vsub5(vdouble v0, vdouble v1, vdouble v2, vdouble v3, vdouble v4) {
  return vsub4(vsub(v0, v1), v2, v3, v4);
}

//

static INLINE vdouble2 normalize_d(vdouble2 t) {
  vdouble2 s;

  s.x = vadd(t.x, t.y);
  s.y = vadd(vsub(t.x, s.x), t.y);

  return s;
}

static INLINE vdouble2 scale_d(vdouble2 d, vdouble s) {
  vdouble2 r = {vmul(d.x, s), vmul(d.y, s)};
  return r;
}

static INLINE vdouble2 add_ss(vdouble x, vdouble y) {
  vdouble2 r;

  r.x = vadd(x, y);
  r.y = vadd(vsub(x, r.x), y);

  return r;
}

static INLINE vdouble2 add2_ss(vdouble x, vdouble y) {
  vdouble2 r;

  r.x = vadd(x, y);
  vdouble v = vsub(r.x, x);
  r.y = vadd(vsub(x, vsub(r.x, v)), vsub(y, v));

  return r;
}

static INLINE vdouble2 add_ds(vdouble2 x, vdouble y) {
  vdouble2 r;

  r.x = vadd(x.x, y);
  r.y = vadd3(vsub(x.x, r.x), y, x.y);

  return r;
}

static INLINE vdouble2 add2_ds(vdouble2 x, vdouble y) {
  vdouble2 r;

  r.x = vadd(x.x, y);
  vdouble v = vsub(r.x, x.x);
  r.y = vadd(vsub(x.x, vsub(r.x, v)), vsub(y, v));
  r.y = vadd(r.y, x.y);

  return r;
}

static INLINE vdouble2 add_sd(vdouble x, vdouble2 y) {
  vdouble2 r;

  r.x = vadd(x, y.x);
  r.y = vadd3(vsub(x, r.x), y.x, y.y);

  return r;
}

static INLINE vdouble2 add_dd(vdouble2 x, vdouble2 y) {
  // |x| >= |y|

  vdouble2 r;

  r.x = vadd(x.x, y.x);
  r.y = vadd4(vsub(x.x, r.x), y.x, x.y, y.y);

  return r;
}

static INLINE vdouble2 add2_dd(vdouble2 x, vdouble2 y) {
  vdouble2 r;

  r.x  = vadd(x.x, y.x);
  vdouble v = vsub(r.x, x.x);
  r.y = vadd(vsub(x.x, vsub(r.x, v)), vsub(y.x, v));
  r.y = vadd(r.y, vadd(x.y, y.y));

  return r;
}

static INLINE vdouble2 div_dd(vdouble2 n, vdouble2 d) {
  vdouble t = vrec(d.x);
  vdouble dh  = vupper(d.x), dl  = vsub(d.x,  dh);
  vdouble th  = vupper(t  ), tl  = vsub(t  ,  th);
  vdouble nhh = vupper(n.x), nhl = vsub(n.x, nhh);

  vdouble2 q;

  q.x = vmul(n.x, t);

  vdouble u = vadd5(vsub(vmul(nhh, th), q.x), vmul(nhh, tl), vmul(nhl, th), vmul(nhl, tl),
		    vmul(q.x, vsub5(vcast_vd_d(1), vmul(dh, th), vmul(dh, tl), vmul(dl, th), vmul(dl, tl))));

  q.y = vadd(vmul(t, vsub(n.y, vmul(q.x, d.y))), u);

  return q;
}

static INLINE vdouble2 mul_ss(vdouble x, vdouble y) {
  vdouble xh = vupper(x), xl = vsub(x, xh);
  vdouble yh = vupper(y), yl = vsub(y, yh);
  vdouble2 r;

  r.x = vmul(x, y);
  r.y = vadd5(vmul(xh, yh), vneg(r.x), vmul(xl, yh), vmul(xh, yl), vmul(xl, yl));

  return r;
}

static INLINE vdouble2 mul_ds(vdouble2 x, vdouble y) {
  vdouble xh = vupper(x.x), xl = vsub(x.x, xh);
  vdouble yh = vupper(y  ), yl = vsub(y  , yh);
  vdouble2 r;

  r.x = vmul(x.x, y);
  r.y = vadd6(vmul(xh, yh), vneg(r.x), vmul(xl, yh), vmul(xh, yl), vmul(xl, yl), vmul(x.y, y));

  return r;
}

static INLINE vdouble2 mul_dd(vdouble2 x, vdouble2 y) {
  vdouble xh = vupper(x.x), xl = vsub(x.x, xh);
  vdouble yh = vupper(y.x), yl = vsub(y.x, yh);
  vdouble2 r;

  r.x = vmul(x.x, y.x);
  r.y = vadd7(vmul(xh, yh), vneg(r.x), vmul(xl, yh), vmul(xh, yl), vmul(xl, yl), vmul(x.x, y.y), vmul(x.y, y.x));

  return r;
}

static INLINE vdouble2 squ_d(vdouble2 x) {
  vdouble xh = vupper(x.x), xl = vsub(x.x, xh);
  vdouble2 r;

  r.x = vmul(x.x, x.x);
  r.y = vadd5(vmul(xh, xh), vneg(r.x), vmul(vadd(xh, xh), xl), vmul(xl, xl), vmul(x.x, vadd(x.y, x.y)));

  return r;
}

static INLINE vdouble2 rec_s(vdouble d) {
  vdouble t = vrec(d);
  vdouble dh = vupper(d), dl = vsub(d, dh);
  vdouble th = vupper(t), tl = vsub(t, th);
  vdouble2 q;

  q.x = t;
  q.y = vmul(t, vsub5(vcast_vd_d(1), vmul(dh, th), vmul(dh, tl), vmul(dl, th), vmul(dl, tl)));

  return q;
}

static INLINE vdouble2 sqrt_d(vdouble2 d) {
  vdouble t = vsqrt(vadd(d.x, d.y));
  return scale_d(mul_dd(add2_dd(d, mul_ss(t, t)), rec_s(t)), vcast_vd_d(0.5));
}

//

static INLINE vdouble xldexp(vdouble x, vint q) { return vldexp(x, q); }

static INLINE vint xilogb(vdouble d) {
  vdouble e = vcast_vd_vi(vsubi(vilogbp1(vabs(d)), vcast_vi_i(1)));
  e = vsel(vmask_eq(d, vcast_vd_d(0)), vcast_vd_d(-2147483648.0), e);
  e = vsel(vmask_eq(vabs(d), vcast_vd_d(rtengine::RT_INFINITY)), vcast_vd_d(2147483647), e);
  return vrint_vi_vd(e);
}

static INLINE vdouble xsin(vdouble d) {
  vint q;
  vdouble u, s;

  q = vrint_vi_vd(vmul(d, vcast_vd_d(rtengine::RT_1_PI)));

  u = vcast_vd_vi(q);
  d = vadd(d, vmul(u, vcast_vd_d(-PI4_A*4)));
  d = vadd(d, vmul(u, vcast_vd_d(-PI4_B*4)));
  d = vadd(d, vmul(u, vcast_vd_d(-PI4_C*4)));

  s = vmul(d, d);

  d = vsel(vmaski_eq(vandi(q, vcast_vi_i(1)), vcast_vi_i(1)), vneg(d), d);

  u = vcast_vd_d(-7.97255955009037868891952e-18);
  u = vmla(u, s, vcast_vd_d(2.81009972710863200091251e-15));
  u = vmla(u, s, vcast_vd_d(-7.64712219118158833288484e-13));
  u = vmla(u, s, vcast_vd_d(1.60590430605664501629054e-10));
  u = vmla(u, s, vcast_vd_d(-2.50521083763502045810755e-08));
  u = vmla(u, s, vcast_vd_d(2.75573192239198747630416e-06));
  u = vmla(u, s, vcast_vd_d(-0.000198412698412696162806809));
  u = vmla(u, s, vcast_vd_d(0.00833333333333332974823815));
  u = vmla(u, s, vcast_vd_d(-0.166666666666666657414808));

  u = vmla(s, vmul(u, d), d);

  return u;
}

static INLINE vdouble xcos(vdouble d) {
  vint q;
  vdouble u, s;

  q = vrint_vi_vd(vsub(vmul(d, vcast_vd_d(rtengine::RT_1_PI)), vcast_vd_d(0.5)));
  q = vaddi(vaddi(q, q), vcast_vi_i(1));

  u = vcast_vd_vi(q);
  d = vadd(d, vmul(u, vcast_vd_d(-PI4_A*2)));
  d = vadd(d, vmul(u, vcast_vd_d(-PI4_B*2)));
  d = vadd(d, vmul(u, vcast_vd_d(-PI4_C*2)));

  s = vmul(d, d);

  d = vsel(vmaski_eq(vandi(q, vcast_vi_i(2)), vcast_vi_i(0)), vneg(d), d);

  u = vcast_vd_d(-7.97255955009037868891952e-18);
  u = vmla(u, s, vcast_vd_d(2.81009972710863200091251e-15));
  u = vmla(u, s, vcast_vd_d(-7.64712219118158833288484e-13));
  u = vmla(u, s, vcast_vd_d(1.60590430605664501629054e-10));
  u = vmla(u, s, vcast_vd_d(-2.50521083763502045810755e-08));
  u = vmla(u, s, vcast_vd_d(2.75573192239198747630416e-06));
  u = vmla(u, s, vcast_vd_d(-0.000198412698412696162806809));
  u = vmla(u, s, vcast_vd_d(0.00833333333333332974823815));
  u = vmla(u, s, vcast_vd_d(-0.166666666666666657414808));

  u = vmla(s, vmul(u, d), d);

  return u;
}

static INLINE vdouble2 xsincos(vdouble d) {
  vint q;
  vmask m;
  vdouble u, s, t, rx, ry;
  vdouble2 r;

  q = vrint_vi_vd(vmul(d, vcast_vd_d(rtengine::RT_2_PI)));

  s = d;

  u = vcast_vd_vi(q);
  s = vmla(u, vcast_vd_d(-PI4_A*2), s);
  s = vmla(u, vcast_vd_d(-PI4_B*2), s);
  s = vmla(u, vcast_vd_d(-PI4_C*2), s);

  t = s;

  s = vmul(s, s);

  u = vcast_vd_d(1.58938307283228937328511e-10);
  u = vmla(u, s, vcast_vd_d(-2.50506943502539773349318e-08));
  u = vmla(u, s, vcast_vd_d(2.75573131776846360512547e-06));
  u = vmla(u, s, vcast_vd_d(-0.000198412698278911770864914));
  u = vmla(u, s, vcast_vd_d(0.0083333333333191845961746));
  u = vmla(u, s, vcast_vd_d(-0.166666666666666130709393));
  u = vmul(vmul(u, s), t);

  rx = vadd(t, u);

  u = vcast_vd_d(-1.13615350239097429531523e-11);
  u = vmla(u, s, vcast_vd_d(2.08757471207040055479366e-09));
  u = vmla(u, s, vcast_vd_d(-2.75573144028847567498567e-07));
  u = vmla(u, s, vcast_vd_d(2.48015872890001867311915e-05));
  u = vmla(u, s, vcast_vd_d(-0.00138888888888714019282329));
  u = vmla(u, s, vcast_vd_d(0.0416666666666665519592062));
  u = vmla(u, s, vcast_vd_d(-0.5));

  ry = vadd(vcast_vd_d(1), vmul(s, u));

  m = vmaski_eq(vandi(q, vcast_vi_i(1)), vcast_vi_i(0));
  r.x = vsel(m, rx, ry);
  r.y = vsel(m, ry, rx);

  m = vmaski_eq(vandi(q, vcast_vi_i(2)), vcast_vi_i(2));
  r.x = vreinterpret_vd_vm(vxorm(vandm(m, vreinterpret_vm_vd(vcast_vd_d(-0.0))), vreinterpret_vm_vd(r.x)));

  m = vmaski_eq(vandi(vaddi(q, vcast_vi_i(1)), vcast_vi_i(2)), vcast_vi_i(2));
  r.y = vreinterpret_vd_vm(vxorm(vandm(m, vreinterpret_vm_vd(vcast_vd_d(-0.0))), vreinterpret_vm_vd(r.y)));

  m = vmask_isinf(d);
  r.x = vsel(m, vcast_vd_d(rtengine::RT_NAN), r.x);
  r.y = vsel(m, vcast_vd_d(rtengine::RT_NAN), r.y);

  return r;
}

static INLINE vdouble xtan(vdouble d) {
  vint q;
  vdouble u, s, x;
  vmask m;

  q = vrint_vi_vd(vmul(d, vcast_vd_d(rtengine::RT_2_PI)));

  u = vcast_vd_vi(q);
  x = vadd(d, vmul(u, vcast_vd_d(-PI4_A*2)));
  x = vadd(x, vmul(u, vcast_vd_d(-PI4_B*2)));
  x = vadd(x, vmul(u, vcast_vd_d(-PI4_C*2)));

  s = vmul(x, x);

  m = vmaski_eq(vandi(q, vcast_vi_i(1)), vcast_vi_i(1));
  x = vsel(m, vneg(x), x);

  u = vcast_vd_d(1.01419718511083373224408e-05);
  u = vmla(u, s, vcast_vd_d(-2.59519791585924697698614e-05));
  u = vmla(u, s, vcast_vd_d(5.23388081915899855325186e-05));
  u = vmla(u, s, vcast_vd_d(-3.05033014433946488225616e-05));
  u = vmla(u, s, vcast_vd_d(7.14707504084242744267497e-05));
  u = vmla(u, s, vcast_vd_d(8.09674518280159187045078e-05));
  u = vmla(u, s, vcast_vd_d(0.000244884931879331847054404));
  u = vmla(u, s, vcast_vd_d(0.000588505168743587154904506));
  u = vmla(u, s, vcast_vd_d(0.00145612788922812427978848));
  u = vmla(u, s, vcast_vd_d(0.00359208743836906619142924));
  u = vmla(u, s, vcast_vd_d(0.00886323944362401618113356));
  u = vmla(u, s, vcast_vd_d(0.0218694882853846389592078));
  u = vmla(u, s, vcast_vd_d(0.0539682539781298417636002));
  u = vmla(u, s, vcast_vd_d(0.133333333333125941821962));
  u = vmla(u, s, vcast_vd_d(0.333333333333334980164153));

  u = vmla(s, vmul(u, x), x);

  u = vsel(m, vrec(u), u);

  u = vsel(vmask_isinf(d), vcast_vd_d(rtengine::RT_NAN), u);

  return u;
}

static INLINE vdouble atan2k(vdouble y, vdouble x) {
  vdouble s, t, u;
  vint q;
  vmask p;

  q = vseli_lt(x, vcast_vd_d(0), vcast_vi_i(-2), vcast_vi_i(0));
  x = vabs(x);

  q = vseli_lt(x, y, vaddi(q, vcast_vi_i(1)), q);
  p = vmask_lt(x, y);
  s = vsel (p, vneg(x), y);
  t = vmax (x, y);

  s = vdiv(s, t);
  t = vmul(s, s);

  u = vcast_vd_d(-1.88796008463073496563746e-05);
  u = vmla(u, t, vcast_vd_d(0.000209850076645816976906797));
  u = vmla(u, t, vcast_vd_d(-0.00110611831486672482563471));
  u = vmla(u, t, vcast_vd_d(0.00370026744188713119232403));
  u = vmla(u, t, vcast_vd_d(-0.00889896195887655491740809));
  u = vmla(u, t, vcast_vd_d(0.016599329773529201970117));
  u = vmla(u, t, vcast_vd_d(-0.0254517624932312641616861));
  u = vmla(u, t, vcast_vd_d(0.0337852580001353069993897));
  u = vmla(u, t, vcast_vd_d(-0.0407629191276836500001934));
  u = vmla(u, t, vcast_vd_d(0.0466667150077840625632675));
  u = vmla(u, t, vcast_vd_d(-0.0523674852303482457616113));
  u = vmla(u, t, vcast_vd_d(0.0587666392926673580854313));
  u = vmla(u, t, vcast_vd_d(-0.0666573579361080525984562));
  u = vmla(u, t, vcast_vd_d(0.0769219538311769618355029));
  u = vmla(u, t, vcast_vd_d(-0.090908995008245008229153));
  u = vmla(u, t, vcast_vd_d(0.111111105648261418443745));
  u = vmla(u, t, vcast_vd_d(-0.14285714266771329383765));
  u = vmla(u, t, vcast_vd_d(0.199999999996591265594148));
  u = vmla(u, t, vcast_vd_d(-0.333333333333311110369124));

  t = vadd(s, vmul(s, vmul(t, u)));
  t = vadd(t, vmul(vcast_vd_vi(q), vcast_vd_d(rtengine::RT_PI/2)));

  return t;
}

static INLINE vdouble xatan2(vdouble y, vdouble x) {
  vdouble r = atan2k(vabs(y), x);

  r = vmulsign(r, x);
  r = vsel(vorm(vmask_isinf(x), vmask_eq(x, vcast_vd_d(0))), vsub(vcast_vd_d(rtengine::RT_PI/2), visinf2(x, vmulsign(vcast_vd_d(rtengine::RT_PI/2), x))), r);
  r = vsel(vmask_isinf(y), vsub(vcast_vd_d(rtengine::RT_PI/2), visinf2(x, vmulsign(vcast_vd_d(rtengine::RT_PI/4), x))), r);
  r = vsel(vmask_eq(y, vcast_vd_d(0)), vsel(vmask_eq(vsign(x), vcast_vd_d(-1.0)), vcast_vd_d(rtengine::RT_PI), vcast_vd_d(0)), r);

  return vsel(vorm(vmask_isnan(x), vmask_isnan(y)), vcast_vd_d(rtengine::RT_NAN), vmulsign(r, y));
}

static INLINE vdouble xasin(vdouble d) {
  vdouble x, y;
  x = vadd(vcast_vd_d(1), d);
  y = vsub(vcast_vd_d(1), d);
  x = vmul(x, y);
  x = vsqrt(x);
  x = vsel(vmask_isnan(x), vcast_vd_d(rtengine::RT_NAN), atan2k(vabs(d), x));
  return vmulsign(x, d);
}

static INLINE vdouble xacos(vdouble d) {
  vdouble x, y;
  x = vadd(vcast_vd_d(1), d);
  y = vsub(vcast_vd_d(1), d);
  x = vmul(x, y);
  x = vsqrt(x);
  x = vmulsign(atan2k(x, vabs(d)), d);
  y = (vdouble)vandm(vmask_lt(d, vcast_vd_d(0)), (vmask)vcast_vd_d(rtengine::RT_PI));
  x = vadd(x, y);
  return x;
}

static INLINE vdouble xatan(vdouble s) {
  vdouble t, u;
  vint q;

  q = vseli_lt(s, vcast_vd_d(0), vcast_vi_i(2), vcast_vi_i(0));
  s = vabs(s);

  q = vseli_lt(vcast_vd_d(1), s, vaddi(q, vcast_vi_i(1)), q);
  s = vsel(vmask_lt(vcast_vd_d(1), s), vdiv(vcast_vd_d(1), s), s);

  t = vmul(s, s);

  u = vcast_vd_d(-1.88796008463073496563746e-05);
  u = vmla(u, t, vcast_vd_d(0.000209850076645816976906797));
  u = vmla(u, t, vcast_vd_d(-0.00110611831486672482563471));
  u = vmla(u, t, vcast_vd_d(0.00370026744188713119232403));
  u = vmla(u, t, vcast_vd_d(-0.00889896195887655491740809));
  u = vmla(u, t, vcast_vd_d(0.016599329773529201970117));
  u = vmla(u, t, vcast_vd_d(-0.0254517624932312641616861));
  u = vmla(u, t, vcast_vd_d(0.0337852580001353069993897));
  u = vmla(u, t, vcast_vd_d(-0.0407629191276836500001934));
  u = vmla(u, t, vcast_vd_d(0.0466667150077840625632675));
  u = vmla(u, t, vcast_vd_d(-0.0523674852303482457616113));
  u = vmla(u, t, vcast_vd_d(0.0587666392926673580854313));
  u = vmla(u, t, vcast_vd_d(-0.0666573579361080525984562));
  u = vmla(u, t, vcast_vd_d(0.0769219538311769618355029));
  u = vmla(u, t, vcast_vd_d(-0.090908995008245008229153));
  u = vmla(u, t, vcast_vd_d(0.111111105648261418443745));
  u = vmla(u, t, vcast_vd_d(-0.14285714266771329383765));
  u = vmla(u, t, vcast_vd_d(0.199999999996591265594148));
  u = vmla(u, t, vcast_vd_d(-0.333333333333311110369124));

  t = vadd(s, vmul(s, vmul(t, u)));

  t = vsel(vmaski_eq(vandi(q, vcast_vi_i(1)), vcast_vi_i(1)), vsub(vcast_vd_d(rtengine::RT_PI/2), t), t);
  t = vsel(vmaski_eq(vandi(q, vcast_vi_i(2)), vcast_vi_i(2)), vneg(t), t);

  return t;
}

static INLINE vdouble xlog(vdouble d) {
  vdouble x, x2;
  vdouble t, m;
  vint e;

  e = vilogbp1(vmul(d, vcast_vd_d(0.7071)));
  m = vldexp(d, vsubi(vcast_vi_i(0), e));

  x = vdiv(vadd(vcast_vd_d(-1), m), vadd(vcast_vd_d(1), m));
  x2 = vmul(x, x);

  t = vcast_vd_d(0.148197055177935105296783);
  t = vmla(t, x2, vcast_vd_d(0.153108178020442575739679));
  t = vmla(t, x2, vcast_vd_d(0.181837339521549679055568));
  t = vmla(t, x2, vcast_vd_d(0.22222194152736701733275));
  t = vmla(t, x2, vcast_vd_d(0.285714288030134544449368));
  t = vmla(t, x2, vcast_vd_d(0.399999999989941956712869));
  t = vmla(t, x2, vcast_vd_d(0.666666666666685503450651));
  t = vmla(t, x2, vcast_vd_d(2));

  x = vadd(vmul(x, t), vmul(vcast_vd_d(0.693147180559945286226764), vcast_vd_vi(e)));

  x = vsel(vmask_ispinf(d), vcast_vd_d(rtengine::RT_INFINITY), x);
  x = vsel(vmask_gt(vcast_vd_d(0), d), vcast_vd_d(rtengine::RT_NAN), x);
  x = vsel(vmask_eq(d, vcast_vd_d(0)), vcast_vd_d(-rtengine::RT_INFINITY), x);

  return x;
}

static INLINE vdouble xexp(vdouble d) {
  vint q = vrint_vi_vd(vmul(d, vcast_vd_d(R_LN2)));
  vdouble s, u;

  s = vadd(d, vmul(vcast_vd_vi(q), vcast_vd_d(-L2U)));
  s = vadd(s, vmul(vcast_vd_vi(q), vcast_vd_d(-L2L)));

  u = vcast_vd_d(2.08860621107283687536341e-09);
  u = vmla(u, s, vcast_vd_d(2.51112930892876518610661e-08));
  u = vmla(u, s, vcast_vd_d(2.75573911234900471893338e-07));
  u = vmla(u, s, vcast_vd_d(2.75572362911928827629423e-06));
  u = vmla(u, s, vcast_vd_d(2.4801587159235472998791e-05));
  u = vmla(u, s, vcast_vd_d(0.000198412698960509205564975));
  u = vmla(u, s, vcast_vd_d(0.00138888888889774492207962));
  u = vmla(u, s, vcast_vd_d(0.00833333333331652721664984));
  u = vmla(u, s, vcast_vd_d(0.0416666666666665047591422));
  u = vmla(u, s, vcast_vd_d(0.166666666666666851703837));
  u = vmla(u, s, vcast_vd_d(0.5));

  u = vadd(vcast_vd_d(1), vadd(s, vmul(vmul(s, s), u)));

  u = vldexp(u, q);

  u = vsel(vmask_isminf(d), vcast_vd_d(0), u);

  return u;
}

static INLINE vdouble2 logk(vdouble d) {
  vdouble2 x, x2;
  vdouble t, m;
  vint e;

  e = vilogbp1(vmul(d, vcast_vd_d(0.7071)));
  m = vldexp(d, vsubi(vcast_vi_i(0), e));

  x = div_dd(add2_ss(vcast_vd_d(-1), m), add2_ss(vcast_vd_d(1), m));
  x2 = squ_d(x);
  x2 = normalize_d(x2);

  t = vcast_vd_d(0.134601987501262130076155);
  t = vmla(t, x2.x, vcast_vd_d(0.132248509032032670243288));
  t = vmla(t, x2.x, vcast_vd_d(0.153883458318096079652524));
  t = vmla(t, x2.x, vcast_vd_d(0.181817427573705403298686));
  t = vmla(t, x2.x, vcast_vd_d(0.222222231326187414840781));
  t = vmla(t, x2.x, vcast_vd_d(0.285714285651261412873718));
  t = vmla(t, x2.x, vcast_vd_d(0.400000000000222439910458));
  t = vmla(t, x2.x, vcast_vd_d(0.666666666666666371239645));

  return add2_dd(mul_ds(dd(vcast_vd_d(0.693147180559945286226764), vcast_vd_d(2.319046813846299558417771e-17)),
		       vcast_vd_vi(e)),
		add2_dd(scale_d(x, vcast_vd_d(2)), mul_ds(mul_dd(x2, x), t)));
}

static INLINE vdouble expk(vdouble2 d) {
  vdouble u = vmul(vadd(d.x, d.y), vcast_vd_d(R_LN2));
  vint q = vrint_vi_vd(u);
  vdouble2 s, t;

  s = add2_ds(d, vmul(vcast_vd_vi(q), vcast_vd_d(-L2U)));
  s = add2_ds(s, vmul(vcast_vd_vi(q), vcast_vd_d(-L2L)));

  q = vrint_vi_vd(vmin(vmax(vcast_vd_d(-2047.49), u), vcast_vd_d(2047.49)));

  s = normalize_d(s);

  u = vcast_vd_d(2.51069683420950419527139e-08);
  u = vmla(u, s.x, vcast_vd_d(2.76286166770270649116855e-07));
  u = vmla(u, s.x, vcast_vd_d(2.75572496725023574143864e-06));
  u = vmla(u, s.x, vcast_vd_d(2.48014973989819794114153e-05));
  u = vmla(u, s.x, vcast_vd_d(0.000198412698809069797676111));
  u = vmla(u, s.x, vcast_vd_d(0.0013888888939977128960529));
  u = vmla(u, s.x, vcast_vd_d(0.00833333333332371417601081));
  u = vmla(u, s.x, vcast_vd_d(0.0416666666665409524128449));
  u = vmla(u, s.x, vcast_vd_d(0.166666666666666740681535));
  u = vmla(u, s.x, vcast_vd_d(0.500000000000000999200722));

  t = add_dd(s, mul_ds(squ_d(s), u));

  t = add_sd(vcast_vd_d(1), t);
  u = vadd(t.x, t.y);
  u = vldexp(u, q);

  return u;
}

static INLINE vdouble xpow(vdouble x, vdouble y) {
#if 1
  vmask yisint = vmask_eq(vcast_vd_vi(vrint_vi_vd(y)), y);
  vmask yisodd = vandm(vmaski_eq(vandi(vrint_vi_vd(y), vcast_vi_i(1)), vcast_vi_i(1)), yisint);

  vdouble result = expk(mul_ds(logk(vabs(x)), y));

  //result = vsel(vmask_isnan(result), vcast_vd_d(rtengine::RT_INFINITY), result);

  result = vmul(result,
		vsel(vmask_gt(x, vcast_vd_d(0)),
		     vcast_vd_d(1),
		     vsel(yisint,
			  vsel(yisodd,
			       vcast_vd_d(-1),
			       vcast_vd_d(1)),
			  vcast_vd_d(rtengine::RT_NAN))));

  vdouble efx = vreinterpret_vd_vm(vxorm(vreinterpret_vm_vd(vsub(vabs(x), vcast_vd_d(1))), vsignbit(y)));

  result = vsel(vmask_isinf(y),
		vsel(vmask_lt(efx, vcast_vd_d(0)),
		     vcast_vd_d(0),
		     vsel(vmask_eq(efx, vcast_vd_d(0)),
			  vcast_vd_d(1.0),
			  vcast_vd_d(rtengine::RT_INFINITY))),
		result);

  result = vsel(vorm(vmask_isinf(x), vmask_eq(x, vcast_vd_d(0))),
		vmul(vsel(yisodd, vsign(x), vcast_vd_d(1)),
		     vsel(vmask_lt(vsel(vmask_eq(x, vcast_vd_d(0)), vneg(y), y), vcast_vd_d(0)),
			  vcast_vd_d(0),
			  vcast_vd_d(rtengine::RT_INFINITY))),
		result);

  result = vsel(vorm(vmask_isnan(x), vmask_isnan(y)), vcast_vd_d(rtengine::RT_NAN), result);

  result = vsel(vorm(vmask_eq(y, vcast_vd_d(0)), vmask_eq(x, vcast_vd_d(1))), vcast_vd_d(1), result);

  return result;
#else
  return expk(mul_ds(logk(x), y));
#endif
}

static INLINE vdouble2 expk2(vdouble2 d) {
  vdouble u = vmul(vadd(d.x, d.y), vcast_vd_d(R_LN2));
  vint q = vrint_vi_vd(u);
  vdouble2 s, t;

  s = add2_ds(d, vmul(vcast_vd_vi(q), vcast_vd_d(-L2U)));
  s = add2_ds(s, vmul(vcast_vd_vi(q), vcast_vd_d(-L2L)));

  q = vrint_vi_vd(vmin(vmax(vcast_vd_d(-2047.49), u), vcast_vd_d(2047.49)));

  s = normalize_d(s);

  u = vcast_vd_d(2.51069683420950419527139e-08);
  u = vmla(u, s.x, vcast_vd_d(2.76286166770270649116855e-07));
  u = vmla(u, s.x, vcast_vd_d(2.75572496725023574143864e-06));
  u = vmla(u, s.x, vcast_vd_d(2.48014973989819794114153e-05));
  u = vmla(u, s.x, vcast_vd_d(0.000198412698809069797676111));
  u = vmla(u, s.x, vcast_vd_d(0.0013888888939977128960529));
  u = vmla(u, s.x, vcast_vd_d(0.00833333333332371417601081));
  u = vmla(u, s.x, vcast_vd_d(0.0416666666665409524128449));
  u = vmla(u, s.x, vcast_vd_d(0.166666666666666740681535));
  u = vmla(u, s.x, vcast_vd_d(0.500000000000000999200722));

  t = add_dd(s, mul_ds(squ_d(s), u));

  t = add_sd(vcast_vd_d(1), t);

  return dd(vldexp(t.x, q), vldexp(t.y, q));
}

static INLINE vdouble xsinh(vdouble x) {
  vdouble y = vabs(x);
  vdouble2 d = expk2(dd(y, vcast_vd_d(0)));
  d = add2_dd(d, div_dd(dd(vcast_vd_d(-1), vcast_vd_d(0)), d));
  y = vmul(vadd(d.x, d.y), vcast_vd_d(0.5));

  y = vsel(vorm(vmask_isinf(x), vmask_isnan(y)), vcast_vd_d(rtengine::RT_INFINITY), y);
  y = vmulsign(y, x);
  y = vsel(vmask_isnan(x), vcast_vd_d(rtengine::RT_NAN), y);

  return y;
}

static INLINE vdouble xcosh(vdouble x) {
  vdouble2 d = expk2(dd(x, vcast_vd_d(0)));
  d = add2_dd(d, div_dd(dd(vcast_vd_d(1), vcast_vd_d(0)), d));
  vdouble y = vmul(vadd(d.x, d.y), vcast_vd_d(0.5));

  y = vsel(vorm(vmask_isinf(x), vmask_isnan(y)), vcast_vd_d(rtengine::RT_INFINITY), y);
  y = vsel(vmask_isnan(x), vcast_vd_d(rtengine::RT_NAN), y);

  return y;
}

static INLINE vdouble xtanh(vdouble x) {
  vdouble y = vabs(x);
  vdouble2 d = expk2(dd(y, vcast_vd_d(0)));
  vdouble2 e = div_dd(dd(vcast_vd_d(1), vcast_vd_d(0)), d);
  d = div_dd(add2_dd(d, scale_d(e, vcast_vd_d(-1))), add2_dd(d, e));
  y = d.x + d.y;

  y = vsel(vorm(vmask_isinf(x), vmask_isnan(y)), vcast_vd_d(1.0), y);
  y = vmulsign(y, x);
  y = vsel(vmask_isnan(x), vcast_vd_d(rtengine::RT_NAN), y);

  return y;
}

static INLINE vdouble2 logk2(vdouble2 d) {
  vdouble2 x, x2, m;
  vdouble t;
  vint e;

  d = normalize_d(d);
  e = vilogbp1(vmul(d.x, vcast_vd_d(0.7071)));
  m = scale_d(d, vldexp(vcast_vd_d(1), vsubi(vcast_vi_i(0), e)));

  x = div_dd(add2_ds(m, vcast_vd_d(-1)), add2_ds(m, vcast_vd_d(1)));
  x2 = squ_d(x);
  x2 = normalize_d(x2);

  t = vcast_vd_d(0.134601987501262130076155);
  t = vmla(t, x2.x, vcast_vd_d(0.132248509032032670243288));
  t = vmla(t, x2.x, vcast_vd_d(0.153883458318096079652524));
  t = vmla(t, x2.x, vcast_vd_d(0.181817427573705403298686));
  t = vmla(t, x2.x, vcast_vd_d(0.222222231326187414840781));
  t = vmla(t, x2.x, vcast_vd_d(0.285714285651261412873718));
  t = vmla(t, x2.x, vcast_vd_d(0.400000000000222439910458));
  t = vmla(t, x2.x, vcast_vd_d(0.666666666666666371239645));

  return add2_dd(mul_ds(dd(vcast_vd_d(0.693147180559945286226764), vcast_vd_d(2.319046813846299558417771e-17)),
		       vcast_vd_vi(e)),
		add2_dd(scale_d(x, vcast_vd_d(2)), mul_ds(mul_dd(x2, x), t)));
}

static INLINE vdouble xasinh(vdouble x) {
  vdouble y = vabs(x);
  vdouble2 d = logk2(add2_ds(sqrt_d(add2_ds(mul_ss(y, y),  vcast_vd_d(1))), y));
  y = vadd(d.x, d.y);

  y = vsel(vorm(vmask_isinf(x), vmask_isnan(y)), vcast_vd_d(rtengine::RT_INFINITY), y);
  y = vmulsign(y, x);
  y = vsel(vmask_isnan(x), vcast_vd_d(rtengine::RT_NAN), y);

  return y;
}

static INLINE vdouble xacosh(vdouble x) {
  vdouble2 d = logk2(add2_ds(sqrt_d(add2_ds(mul_ss(x, x), vcast_vd_d(-1))), x));
  vdouble y = vadd(d.x, d.y);

  y = vsel(vorm(vmask_isinf(x), vmask_isnan(y)), vcast_vd_d(rtengine::RT_INFINITY), y);
  y = vsel(vmask_eq(x, vcast_vd_d(1.0)), vcast_vd_d(0.0), y);
  y = vsel(vmask_lt(x, vcast_vd_d(1.0)), vcast_vd_d(rtengine::RT_NAN), y);
  y = vsel(vmask_isnan(x), vcast_vd_d(rtengine::RT_NAN), y);

  return y;
}

static INLINE vdouble xatanh(vdouble x) {
  vdouble y = vabs(x);
  vdouble2 d = logk2(div_dd(add2_ss(vcast_vd_d(1), y), add2_ss(vcast_vd_d(1), -y)));
  y = vsel(vmask_gt(y, vcast_vd_d(1.0)), vcast_vd_d(rtengine::RT_NAN), vsel(vmask_eq(y, vcast_vd_d(1.0)), vcast_vd_d(rtengine::RT_INFINITY), vmul(vadd(d.x, d.y), vcast_vd_d(0.5))));

  y = vsel(vorm(vmask_isinf(x), vmask_isnan(y)), vcast_vd_d(rtengine::RT_NAN), y);
  y = vmulsign(y, x);
  y = vsel(vmask_isnan(x), vcast_vd_d(rtengine::RT_NAN), y);

  return y;
}

static INLINE vdouble xcbrt(vdouble d) {
  vdouble x, y, q = vcast_vd_d(1.0);
  vint e, qu, re;
  vdouble t;

  e = vilogbp1(vabs(d));
  d = vldexp(d, vsubi(vcast_vi_i(0), e));

  t = vadd(vcast_vd_vi(e), vcast_vd_d(6144));
  qu = vtruncate_vi_vd(vdiv(t, vcast_vd_d(3)));
  re = vtruncate_vi_vd(vsub(t, vmul(vcast_vd_vi(qu), vcast_vd_d(3))));

  q = vsel(vmaski_eq(re, vcast_vi_i(1)), vcast_vd_d(1.2599210498948731647672106), q);
  q = vsel(vmaski_eq(re, vcast_vi_i(2)), vcast_vd_d(1.5874010519681994747517056), q);
  q = vldexp(q, vsubi(qu, vcast_vi_i(2048)));

  q = vmulsign(q, d);

  d = vabs(d);

  x = vcast_vd_d(-0.640245898480692909870982);
  x = vmla(x, d, vcast_vd_d(2.96155103020039511818595));
  x = vmla(x, d, vcast_vd_d(-5.73353060922947843636166));
  x = vmla(x, d, vcast_vd_d(6.03990368989458747961407));
  x = vmla(x, d, vcast_vd_d(-3.85841935510444988821632));
  x = vmla(x, d, vcast_vd_d(2.2307275302496609725722));

  y = vmul(x, x); y = vmul(y, y); x = vsub(x, vmul(vmla(d, y, vneg(x)), vcast_vd_d(1.0 / 3.0)));
  y = vmul(vmul(d, x), x);
  y = vmul(vsub(y, vmul(vmul(vcast_vd_d(2.0 / 3.0), y), vmla(y, x, vcast_vd_d(-1.0)))), q);

  return y;
}

static INLINE vdouble xexp2(vdouble a) {
  vdouble u = expk(mul_ds(dd(vcast_vd_d(0.69314718055994528623), vcast_vd_d(2.3190468138462995584e-17)), a));
  u = vsel(vmask_ispinf(a), vcast_vd_d(rtengine::RT_INFINITY), u);
  u = vsel(vmask_isminf(a), vcast_vd_d(0), u);
  return u;
}

static INLINE vdouble xexp10(vdouble a) {
  vdouble u = expk(mul_ds(dd(vcast_vd_d(2.3025850929940459011), vcast_vd_d(-2.1707562233822493508e-16)), a));
  u = vsel(vmask_ispinf(a), vcast_vd_d(rtengine::RT_INFINITY), u);
  u = vsel(vmask_isminf(a), vcast_vd_d(0), u);
  return u;
}

static INLINE vdouble xexpm1(vdouble a) {
  vdouble2 d = add2_ds(expk2(dd(a, vcast_vd_d(0))), vcast_vd_d(-1.0));
  vdouble x = d.x + d.y;
  x = vsel(vmask_ispinf(a), vcast_vd_d(rtengine::RT_INFINITY), x);
  x = vsel(vmask_isminf(a), vcast_vd_d(-1), x);
  return x;
}

static INLINE vdouble xlog10(vdouble a) {
  vdouble2 d = mul_dd(logk(a), dd(vcast_vd_d(0.43429448190325176116), vcast_vd_d(6.6494347733425473126e-17)));
  vdouble x = d.x + d.y;

  x = vsel(vmask_ispinf(a), vcast_vd_d(rtengine::RT_INFINITY), x);
  x = vsel(vmask_gt(vcast_vd_d(0), a), vcast_vd_d(rtengine::RT_NAN), x);
  x = vsel(vmask_eq(a, vcast_vd_d(0)), vcast_vd_d(-rtengine::RT_INFINITY), x);

  return x;
}

static INLINE vdouble xlog1p(vdouble a) {
  vdouble2 d = logk2(add2_ss(a, vcast_vd_d(1)));
  vdouble x = d.x + d.y;

  x = vsel(vmask_ispinf(a), vcast_vd_d(rtengine::RT_INFINITY), x);
  x = vsel(vmask_gt(vcast_vd_d(-1), a), vcast_vd_d(rtengine::RT_NAN), x);
  x = vsel(vmask_eq(a, vcast_vd_d(-1)), vcast_vd_d(-rtengine::RT_INFINITY), x);

  return x;
}

//

typedef struct {
  vfloat x, y;
} vfloat2;

static INLINE vfloat vabsf(vfloat f) { return (vfloat)vandnotm((vmask)vcast_vf_f(-0.0f), (vmask)f); }
static INLINE vfloat vnegf(vfloat f) { return (vfloat)vxorm((vmask)f, (vmask)vcast_vf_f(-0.0f)); }

#if defined( __SSE4_1__ ) && defined( __x86_64__ )
	// only one instruction when using SSE4.1
	static INLINE vfloat vself(vmask mask, vfloat x, vfloat y) {
		return _mm_blendv_ps(y,x,(vfloat)mask);
	}

	static INLINE vint vselc(vmask mask, vint x, vint y) {
		return _mm_blendv_epi8(y,x,mask);
	}

#else
	// three instructions when using SSE2
	static INLINE vfloat vself(vmask mask, vfloat x, vfloat y) {
		return (vfloat)vorm(vandm(mask, (vmask)x), vandnotm(mask, (vmask)y));
	}

	static INLINE vint vselc(vmask mask, vint x, vint y) {
	    return vorm(vandm(mask, (vmask)x), vandnotm(mask, (vmask)y));
	}
#endif

static INLINE vfloat vselfzero(vmask mask, vfloat x) {
     // returns value of x if corresponding mask bits are 1, else returns 0
     // faster than vself(mask, x, ZEROV)
    return _mm_and_ps((vfloat)mask, x);
}
static INLINE vfloat vselfnotzero(vmask mask, vfloat x) {
    // returns value of x if corresponding mask bits are 0, else returns 0
    // faster than vself(mask, ZEROV, x)
    return _mm_andnot_ps((vfloat)mask, x);
}

static INLINE vint vselizero(vmask mask, vint x) {
     // returns value of x if corresponding mask bits are 1, else returns 0
     // faster than vselc(mask, x, ZEROV)
    return _mm_and_si128(mask, x);
}
static INLINE vint vselinotzero(vmask mask, vint x) {
    // returns value of x if corresponding mask bits are 0, else returns 0
    // faster than vselc(mask, ZEROV, x)
    return _mm_andnot_si128(mask, x);
}

static INLINE vint2 vseli2_lt(vfloat f0, vfloat f1, vint2 x, vint2 y) {
  vint2 m2 = vcast_vi2_vm(vmaskf_lt(f0, f1));
  return vori2(vandi2(m2, x), vandnoti2(m2, y));
}

static INLINE vmask vsignbitf(vfloat f) {
  return vandm((vmask)f, (vmask)vcast_vf_f(-0.0f));
}

static INLINE vfloat vmulsignf(vfloat x, vfloat y) {
  return (vfloat)vxorm((vmask)x, vsignbitf(y));
}

static INLINE vfloat vsignf(vfloat f) {
  return (vfloat)vorm((vmask)vcast_vf_f(1.0f), vandm((vmask)vcast_vf_f(-0.0f), (vmask)f));
}

static INLINE vmask vmaskf_isinf(vfloat d) { return vmaskf_eq(vabsf(d), vcast_vf_f(INFINITYf)); }
static INLINE vmask vmaskf_ispinf(vfloat d) { return vmaskf_eq(d, vcast_vf_f(INFINITYf)); }
static INLINE vmask vmaskf_isminf(vfloat d) { return vmaskf_eq(d, vcast_vf_f(-INFINITYf)); }
static INLINE vmask vmaskf_isnan(vfloat d) { return vmaskf_neq(d, d); }
static INLINE vfloat visinf2f(vfloat d, vfloat m) { return (vfloat)vandm(vmaskf_isinf(d), vorm(vsignbitf(d), (vmask)m)); }
static INLINE vfloat visinff(vfloat d) { return visinf2f(d, vcast_vf_f(1.0f)); }

static INLINE vint2 vilogbp1f(vfloat d) {
  vmask m = vmaskf_lt(d, vcast_vf_f(5.421010862427522E-20f));
  d = vself(m, vmulf(vcast_vf_f(1.8446744073709552E19f), d), d);
  vint2 q = vandi2(vsrli2(vcast_vi2_vm(vreinterpret_vm_vf(d)), 23), vcast_vi2_i(0xff));
  q = vsubi2(q, vseli2(m, vcast_vi2_i(64 + 0x7e), vcast_vi2_i(0x7e)));
  return q;
}

static INLINE vfloat vldexpf(vfloat x, vint2 q) {
  vfloat u;
  vint2 m = vsrai2(q, 31);
  m = vslli2(vsubi2(vsrai2(vaddi2(m, q), 6), m), 4);
  q = vsubi2(q, vslli2(m, 2));
  u = vreinterpret_vf_vm(vcast_vm_vi2(vslli2(vaddi2(m, vcast_vi2_i(0x7f)), 23)));
  x = vmulf(vmulf(vmulf(vmulf(x, u), u), u), u);
  u = vreinterpret_vf_vm(vcast_vm_vi2(vslli2(vaddi2(q, vcast_vi2_i(0x7f)), 23)));
  return vmulf(x, u);
}

static INLINE vfloat xsinf(vfloat d) {
  vint2 q;
  vfloat u, s;

  q = vrint_vi2_vf(vmulf(d, vcast_vf_f((float)rtengine::RT_1_PI)));

  u = vcast_vf_vi2(q);
  d = vmlaf(u, vcast_vf_f(-PI4_Af*4), d);
  d = vmlaf(u, vcast_vf_f(-PI4_Bf*4), d);
  d = vmlaf(u, vcast_vf_f(-PI4_Cf*4), d);
  d = vmlaf(u, vcast_vf_f(-PI4_Df*4), d);

  s = vmulf(d, d);

  d = vself(vmaski2_eq(vandi2(q, vcast_vi2_i(1)), vcast_vi2_i(1)), vnegf(d), d);

  u = vcast_vf_f(2.6083159809786593541503e-06f);
  u = vmlaf(u, s, vcast_vf_f(-0.0001981069071916863322258f));
  u = vmlaf(u, s, vcast_vf_f(0.00833307858556509017944336f));
  u = vmlaf(u, s, vcast_vf_f(-0.166666597127914428710938f));

  u = vmlaf(s, vmulf(u, d), d);

  return u;
}

static INLINE vfloat xcosf(vfloat d) {
  vint2 q;
  vfloat u, s;

  q = vrint_vi2_vf(vsubf(vmulf(d, vcast_vf_f((float)rtengine::RT_1_PI)), vcast_vf_f(0.5f)));
  q = vaddi2(vaddi2(q, q), vcast_vi2_i(1));

  u = vcast_vf_vi2(q);
  d = vmlaf(u, vcast_vf_f(-PI4_Af*2), d);
  d = vmlaf(u, vcast_vf_f(-PI4_Bf*2), d);
  d = vmlaf(u, vcast_vf_f(-PI4_Cf*2), d);
  d = vmlaf(u, vcast_vf_f(-PI4_Df*2), d);

  s = vmulf(d, d);

  d = vself(vmaski2_eq(vandi2(q, vcast_vi2_i(2)), vcast_vi2_i(2)), d, vnegf(d));

  u = vcast_vf_f(2.6083159809786593541503e-06f);
  u = vmlaf(u, s, vcast_vf_f(-0.0001981069071916863322258f));
  u = vmlaf(u, s, vcast_vf_f(0.00833307858556509017944336f));
  u = vmlaf(u, s, vcast_vf_f(-0.166666597127914428710938f));

  u = vmlaf(s, vmulf(u, d), d);

  return u;
}

static INLINE vfloat2 xsincosf(vfloat d) {
  vint2 q;
  vmask m;
  vfloat u, s, t, rx, ry;
  vfloat2 r;

  q = vrint_vi2_vf(vmulf(d, vcast_vf_f((float)rtengine::RT_2_PI)));

  s = d;

  u = vcast_vf_vi2(q);
  s = vmlaf(u, vcast_vf_f(-PI4_Af*2), s);
  s = vmlaf(u, vcast_vf_f(-PI4_Bf*2), s);
  s = vmlaf(u, vcast_vf_f(-PI4_Cf*2), s);
  s = vmlaf(u, vcast_vf_f(-PI4_Df*2), s);

  t = s;

  s = vmulf(s, s);

  u = vcast_vf_f(-0.000195169282960705459117889f);
  u = vmlaf(u, s, vcast_vf_f(0.00833215750753879547119141f));
  u = vmlaf(u, s, vcast_vf_f(-0.166666537523269653320312f));
  u = vmulf(vmulf(u, s), t);

  rx = vaddf(t, u);

  u = vcast_vf_f(-2.71811842367242206819355e-07f);
  u = vmlaf(u, s, vcast_vf_f(2.47990446951007470488548e-05f));
  u = vmlaf(u, s, vcast_vf_f(-0.00138888787478208541870117f));
  u = vmlaf(u, s, vcast_vf_f(0.0416666641831398010253906f));
  u = vmlaf(u, s, vcast_vf_f(-0.5));

  ry = vaddf(vcast_vf_f(1), vmulf(s, u));

  m = vmaski2_eq(vandi2(q, vcast_vi2_i(1)), vcast_vi2_i(0));
  r.x = vself(m, rx, ry);
  r.y = vself(m, ry, rx);

  m = vmaski2_eq(vandi2(q, vcast_vi2_i(2)), vcast_vi2_i(2));
  r.x = vreinterpret_vf_vm(vxorm(vandm(m, vreinterpret_vm_vf(vcast_vf_f(-0.0))), vreinterpret_vm_vf(r.x)));

  m = vmaski2_eq(vandi2(vaddi2(q, vcast_vi2_i(1)), vcast_vi2_i(2)), vcast_vi2_i(2));
  r.y = vreinterpret_vf_vm(vxorm(vandm(m, vreinterpret_vm_vf(vcast_vf_f(-0.0))), vreinterpret_vm_vf(r.y)));

  m = vmaskf_isinf(d);
  r.x = vself(m, vcast_vf_f(rtengine::RT_NAN), r.x);
  r.y = vself(m, vcast_vf_f(rtengine::RT_NAN), r.y);

  return r;
}

static INLINE vfloat xtanf(vfloat d) {
  vint2 q;
  vmask m;
  vfloat u, s, x;

  q = vrint_vi2_vf(vmulf(d, vcast_vf_f((float)(2 * rtengine::RT_1_PI))));

  x = d;

  u = vcast_vf_vi2(q);
  x = vmlaf(u, vcast_vf_f(-PI4_Af*2), x);
  x = vmlaf(u, vcast_vf_f(-PI4_Bf*2), x);
  x = vmlaf(u, vcast_vf_f(-PI4_Cf*2), x);
  x = vmlaf(u, vcast_vf_f(-PI4_Df*2), x);

  s = vmulf(x, x);

  m = vmaski2_eq(vandi2(q, vcast_vi2_i(1)), vcast_vi2_i(1));
  x = vself(m, vnegf(x), x);

  u = vcast_vf_f(0.00927245803177356719970703f);
  u = vmlaf(u, s, vcast_vf_f(0.00331984995864331722259521f));
  u = vmlaf(u, s, vcast_vf_f(0.0242998078465461730957031f));
  u = vmlaf(u, s, vcast_vf_f(0.0534495301544666290283203f));
  u = vmlaf(u, s, vcast_vf_f(0.133383005857467651367188f));
  u = vmlaf(u, s, vcast_vf_f(0.333331853151321411132812f));

  u = vmlaf(s, vmulf(u, x), x);

  u = vself(m, vrecf(u), u);

  u = vself(vmaskf_isinf(d), vcast_vf_f(NANf), u);

  return u;
}

static INLINE vfloat xatanf(vfloat s) {
  vfloat t, u;
  vint2 q;

  q = vseli2_lt(s, vcast_vf_f(0.0f), vcast_vi2_i(2), vcast_vi2_i(0));
  s = vabsf(s);

  q = vseli2_lt(vcast_vf_f(1.0f), s, vaddi2(q, vcast_vi2_i(1)), q);
  s = vself(vmaskf_lt(vcast_vf_f(1.0f), s), vdivf(vcast_vf_f(1.0f), s), s);

  t = vmulf(s, s);

  u = vcast_vf_f(0.00282363896258175373077393f);
  u = vmlaf(u, t, vcast_vf_f(-0.0159569028764963150024414f));
  u = vmlaf(u, t, vcast_vf_f(0.0425049886107444763183594f));
  u = vmlaf(u, t, vcast_vf_f(-0.0748900920152664184570312f));
  u = vmlaf(u, t, vcast_vf_f(0.106347933411598205566406f));
  u = vmlaf(u, t, vcast_vf_f(-0.142027363181114196777344f));
  u = vmlaf(u, t, vcast_vf_f(0.199926957488059997558594f));
  u = vmlaf(u, t, vcast_vf_f(-0.333331018686294555664062f));

  t = vaddf(s, vmulf(s, vmulf(t, u)));

  t = vself(vmaski2_eq(vandi2(q, vcast_vi2_i(1)), vcast_vi2_i(1)), vsubf(vcast_vf_f((float)(rtengine::RT_PI/2)), t), t);
  t = vself(vmaski2_eq(vandi2(q, vcast_vi2_i(2)), vcast_vi2_i(2)), vnegf(t), t);

  return t;
}

static INLINE vfloat atan2kf(vfloat y, vfloat x) {
  vfloat s, t, u;
  vint2 q;
  vmask p;

  q = vseli2_lt(x, vcast_vf_f(0.0f), vcast_vi2_i(-2), vcast_vi2_i(0));
  x = vabsf(x);

  q = vseli2_lt(x, y, vaddi2(q, vcast_vi2_i(1)), q);
  p = vmaskf_lt(x, y);
  s = vself(p, vnegf(x), y);
  t = vmaxf(x, y);

  s = vdivf(s, t);
  t = vmulf(s, s);

  u = vcast_vf_f(0.00282363896258175373077393f);
  u = vmlaf(u, t, vcast_vf_f(-0.0159569028764963150024414f));
  u = vmlaf(u, t, vcast_vf_f(0.0425049886107444763183594f));
  u = vmlaf(u, t, vcast_vf_f(-0.0748900920152664184570312f));
  u = vmlaf(u, t, vcast_vf_f(0.106347933411598205566406f));
  u = vmlaf(u, t, vcast_vf_f(-0.142027363181114196777344f));
  u = vmlaf(u, t, vcast_vf_f(0.199926957488059997558594f));
  u = vmlaf(u, t, vcast_vf_f(-0.333331018686294555664062f));

  t = vaddf(s, vmulf(s, vmulf(t, u)));
  t = vaddf(t, vmulf(vcast_vf_vi2(q), vcast_vf_f((float)(rtengine::RT_PI/2))));

  return t;
}

static INLINE vfloat xatan2f(vfloat y, vfloat x) {
  vfloat r = atan2kf(vabsf(y), x);

  r = vmulsignf(r, x);
  r = vself(vorm(vmaskf_isinf(x), vmaskf_eq(x, vcast_vf_f(0.0f))), vsubf(vcast_vf_f((float)(rtengine::RT_PI/2)), visinf2f(x, vmulsignf(vcast_vf_f((float)(rtengine::RT_PI/2)), x))), r);
  r = vself(vmaskf_isinf(y), vsubf(vcast_vf_f((float)(rtengine::RT_PI/2)), visinf2f(x, vmulsignf(vcast_vf_f((float)(rtengine::RT_PI/4)), x))), r);
  r = vself(vmaskf_eq(y, vcast_vf_f(0.0f)), vselfzero(vmaskf_eq(vsignf(x), vcast_vf_f(-1.0f)), vcast_vf_f((float)rtengine::RT_PI)), r);

  return vself(vorm(vmaskf_isnan(x), vmaskf_isnan(y)), vcast_vf_f(NANf), vmulsignf(r, y));
}

static INLINE vfloat xasinf(vfloat d) {
  vfloat x, y;
  x = vaddf(vcast_vf_f(1.0f), d);
  y = vsubf(vcast_vf_f(1.0f), d);
  x = vmulf(x, y);
  x = vsqrtf(x);
  x = vself(vmaskf_isnan(x), vcast_vf_f(NANf), atan2kf(vabsf(d), x));
  return vmulsignf(x, d);
}

static INLINE vfloat xacosf(vfloat d) {
  vfloat x, y;
  x = vaddf(vcast_vf_f(1.0f), d);
  y = vsubf(vcast_vf_f(1.0f), d);
  x = vmulf(x, y);
  x = vsqrtf(x);
  x = vmulsignf(atan2kf(x, vabsf(d)), d);
  y = (vfloat)vandm(vmaskf_lt(d, vcast_vf_f(0.0f)), (vmask)vcast_vf_f((float)rtengine::RT_PI));
  x = vaddf(x, y);
  return x;
}

static INLINE vfloat xlogf(vfloat d) {
  vfloat x, x2, t, m;
  vint2 e;

  e = vilogbp1f(vmulf(d, vcast_vf_f(0.7071f)));
  m = vldexpf(d, vsubi2(vcast_vi2_i(0), e));

  x = vdivf(vaddf(vcast_vf_f(-1.0f), m), vaddf(vcast_vf_f(1.0f), m));
  x2 = vmulf(x, x);

  t = vcast_vf_f(0.2371599674224853515625f);
  t = vmlaf(t, x2, vcast_vf_f(0.285279005765914916992188f));
  t = vmlaf(t, x2, vcast_vf_f(0.400005519390106201171875f));
  t = vmlaf(t, x2, vcast_vf_f(0.666666567325592041015625f));
  t = vmlaf(t, x2, vcast_vf_f(2.0f));

  x = vaddf(vmulf(x, t), vmulf(vcast_vf_f(0.693147180559945286226764f), vcast_vf_vi2(e)));

  x = vself(vmaskf_ispinf(d), vcast_vf_f(INFINITYf), x);
  x = vself(vmaskf_gt(vcast_vf_f(0), d), vcast_vf_f(NANf), x);
  x = vself(vmaskf_eq(d, vcast_vf_f(0)), vcast_vf_f(-INFINITYf), x);

  return x;
}

static INLINE vfloat xlogf0(vfloat d) {
  vfloat x, x2, t, m;
  vint2 e;

  e = vilogbp1f(vmulf(d, vcast_vf_f(0.7071f)));
  m = vldexpf(d, vsubi2(vcast_vi2_i(0), e));

  x = vdivf(vaddf(vcast_vf_f(-1.0f), m), vaddf(vcast_vf_f(1.0f), m));
  x2 = vmulf(x, x);

  t = vcast_vf_f(0.2371599674224853515625f);
  t = vmlaf(t, x2, vcast_vf_f(0.285279005765914916992188f));
  t = vmlaf(t, x2, vcast_vf_f(0.400005519390106201171875f));
  t = vmlaf(t, x2, vcast_vf_f(0.666666567325592041015625f));
  t = vmlaf(t, x2, vcast_vf_f(2.0f));

  x = vaddf(vmulf(x, t), vmulf(vcast_vf_f(0.693147180559945286226764f), vcast_vf_vi2(e)));

  x = vself(vmaskf_ispinf(d), vcast_vf_f(0), x);
  x = vself(vmaskf_gt(vcast_vf_f(0), d), vcast_vf_f(0), x);
  x = vself(vmaskf_eq(d, vcast_vf_f(0)), vcast_vf_f(0), x);

  return x;
}

static INLINE vfloat xlogfNoCheck(vfloat d) { // this version does not check input values. Use it only when you know the input values are > 0 e.g. when filling a lookup table
  vfloat x, x2, t, m;
  vint2 e;

  e = vilogbp1f(vmulf(d, vcast_vf_f(0.7071f)));
  m = vldexpf(d, vsubi2(vcast_vi2_i(0), e));

  x = vdivf(vaddf(vcast_vf_f(-1.0f), m), vaddf(vcast_vf_f(1.0f), m));
  x2 = vmulf(x, x);

  t = vcast_vf_f(0.2371599674224853515625f);
  t = vmlaf(t, x2, vcast_vf_f(0.285279005765914916992188f));
  t = vmlaf(t, x2, vcast_vf_f(0.400005519390106201171875f));
  t = vmlaf(t, x2, vcast_vf_f(0.666666567325592041015625f));
  t = vmlaf(t, x2, vcast_vf_f(2.0f));

  return vaddf(vmulf(x, t), vmulf(vcast_vf_f(0.693147180559945286226764f), vcast_vf_vi2(e)));

}

static INLINE vfloat xexpf(vfloat d) {
  vint2 q = vrint_vi2_vf(vmulf(d, vcast_vf_f(R_LN2f)));
  vfloat s, u;

  s = vmlaf(vcast_vf_vi2(q), vcast_vf_f(-L2Uf),d);
  s = vmlaf(vcast_vf_vi2(q), vcast_vf_f(-L2Lf),s);

  u = vcast_vf_f(0.00136324646882712841033936f);
  u = vmlaf(u, s, vcast_vf_f(0.00836596917361021041870117f));
  u = vmlaf(u, s, vcast_vf_f(0.0416710823774337768554688f));
  u = vmlaf(u, s, vcast_vf_f(0.166665524244308471679688f));
  u = vmlaf(u, s, vcast_vf_f(0.499999850988388061523438f));

  u = vaddf(vcast_vf_f(1.0f), vmlaf(vmulf(s, s), u, s));

  u = vldexpf(u, q);

  u = vself(vmaskf_isminf(d), vcast_vf_f(0.0f), u);
// -104.0
  u = vself(vmaskf_gt(vcast_vf_f(-104), d), vcast_vf_f(0), u);
  return u;
}

static INLINE vfloat xexpfNoCheck(vfloat d) { // this version does not check input values. Use it only when you know the input values are > -104.f e.g. when filling a lookup table
  vint2 q = vrint_vi2_vf(vmulf(d, vcast_vf_f(R_LN2f)));
  vfloat s, u;

  s = vmlaf(vcast_vf_vi2(q), vcast_vf_f(-L2Uf),d);
  s = vmlaf(vcast_vf_vi2(q), vcast_vf_f(-L2Lf),s);

  u = vcast_vf_f(0.00136324646882712841033936f);
  u = vmlaf(u, s, vcast_vf_f(0.00836596917361021041870117f));
  u = vmlaf(u, s, vcast_vf_f(0.0416710823774337768554688f));
  u = vmlaf(u, s, vcast_vf_f(0.166665524244308471679688f));
  u = vmlaf(u, s, vcast_vf_f(0.499999850988388061523438f));

  u = vaddf(vcast_vf_f(1.0f), vmlaf(vmulf(s, s), u, s));

  return vldexpf(u, q);
}

static INLINE vfloat xcbrtf(vfloat d) {
  vfloat x, y, q = vcast_vf_f(1.0), t;
  vint2 e, qu, re;

  e = vilogbp1f(vabsf(d));
  d = vldexpf(d, vsubi2(vcast_vi2_i(0), e));

  t = vaddf(vcast_vf_vi2(e), vcast_vf_f(6144));
  qu = vtruncate_vi2_vf(vdivf(t, vcast_vf_f(3)));
  re = vtruncate_vi2_vf(vsubf(t, vmulf(vcast_vf_vi2(qu), vcast_vf_f(3))));

  q = vself(vmaski2_eq(re, vcast_vi2_i(1)), vcast_vf_f(1.2599210498948731647672106f), q);
  q = vself(vmaski2_eq(re, vcast_vi2_i(2)), vcast_vf_f(1.5874010519681994747517056f), q);
  q = vldexpf(q, vsubi2(qu, vcast_vi2_i(2048)));

  q = vmulsignf(q, d);
  d = vabsf(d);

  x = vcast_vf_f(-0.601564466953277587890625f);
  x = vmlaf(x, d, vcast_vf_f(2.8208892345428466796875f));
  x = vmlaf(x, d, vcast_vf_f(-5.532182216644287109375f));
  x = vmlaf(x, d, vcast_vf_f(5.898262500762939453125f));
  x = vmlaf(x, d, vcast_vf_f(-3.8095417022705078125f));
  x = vmlaf(x, d, vcast_vf_f(2.2241256237030029296875f));

  y = vmulf(vmulf(d, x), x);
  y = vmulf(vsubf(y, vmulf(vmulf(vcast_vf_f(2.0f / 3.0f), y), vmlaf(y, x, vcast_vf_f(-1.0f)))), q);

  return y;
}

static INLINE vfloat LIMV( vfloat a, vfloat b, vfloat c ) {
return vmaxf( b, vminf(a,c));
}

static INLINE vfloat SQRV(vfloat a){
	return a * a;
}

static inline void vswap( vmask condition, vfloat &a, vfloat &b) {
    // conditional swap the elements of two vfloats
    vfloat temp = vself(condition, a, b); // the values which fit to condition
    condition = vnotm(condition); // invert the condition
    a = vself(condition, a, b); // the values which fit to inverted condition
    b = temp;
}

static inline float vhadd( vfloat a ) {
    // returns a[0] + a[1] + a[2] + a[3]
    a += _mm_movehl_ps(a, a);
    return _mm_cvtss_f32(_mm_add_ss(a, _mm_shuffle_ps(a, a, 1)));
}

static INLINE vfloat vmul2f(vfloat a){
    // fastest way to multiply by 2
	return a + a;
}

static INLINE vfloat vintpf(vfloat a, vfloat b, vfloat c) {
    // calculate a * b + (1 - a) * c (interpolate two values)
    // following is valid:
    // vintpf(a, b+x, c+x) = vintpf(a, b, c) + x
    // vintpf(a, b*x, c*x) = vintpf(a, b, c) * x
    return a * (b-c) + c;
}

static INLINE vfloat vdup(vfloat a){
    // returns { a[0],a[0],a[1],a[1] }
    return _mm_unpacklo_ps( a, a );
}

static INLINE vfloat vaddc2vfu(float &a)
{
    // loads a[0]..a[7] and returns { a[0]+a[1], a[2]+a[3], a[4]+a[5], a[6]+a[7] }
    vfloat a1 = _mm_loadu_ps( &a );
    vfloat a2 = _mm_loadu_ps( (&a) + 4 );
    return _mm_shuffle_ps(a1,a2,_MM_SHUFFLE( 2,0,2,0 )) + _mm_shuffle_ps(a1,a2,_MM_SHUFFLE( 3,1,3,1 ));
}

static INLINE vfloat vadivapb (vfloat a, vfloat b) {
    return a / (a+b);
}

static INLINE void vconvertrgbrgbrgbrgb2rrrrggggbbbb (const float * src, vfloat &rv, vfloat &gv, vfloat &bv) { // cool function name, isn't it ? :P
    // converts a sequence of 4 float RGB triplets to 3 red, green and blue quadruples
    rv = _mm_setr_ps(src[0],src[3],src[6],src[9]);
    gv = _mm_setr_ps(src[1],src[4],src[7],src[10]);
    bv = _mm_setr_ps(src[2],src[5],src[8],src[11]);
}

#endif // __SSE2__
#endif // SLEEFSSEAVX
