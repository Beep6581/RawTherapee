////////////////////////////////////////////////////////////////
//
//  this code was taken from http://shibatch.sourceforge.net/
//  Many thanks to the author of original version: Naoki Shibata
//
//  This version contains modifications made by Ingo Weyrich
//
////////////////////////////////////////////////////////////////

#ifndef _SLEEFC_
#define _SLEEFC_

#include <assert.h>
#include <stdint.h>
//#include <math.h>
#include "rt_math.h"
//#include <bits/nan.h>
//#include <bits/inf.h>

#define PI4_A .7853981554508209228515625
#define PI4_B .794662735614792836713604629039764404296875e-8
#define PI4_C .306161699786838294306516483068750264552437361480769e-16
#define M_4_PI 1.273239544735162542821171882678754627704620361328125

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

__inline int xisnan(double x) { return x != x; }
__inline int xisinf(double x) { return x == rtengine::RT_INFINITY || x == -rtengine::RT_INFINITY; }
__inline int xisminf(double x) { return x == -rtengine::RT_INFINITY; }
__inline int xispinf(double x) { return x == rtengine::RT_INFINITY; }

__inline double ldexpk(double x, int q) {
  double u;
  int m;
  m = q >> 31;
  m = (((m + q) >> 9) - m) << 7;
  q = q - (m << 2);
  u = longBitsToDouble(((int64_t)(m + 0x3ff)) << 52);
  double u2 = u*u;
  u2 = u2 * u2;
  x = x * u2;
  u = longBitsToDouble(((int64_t)(q + 0x3ff)) << 52);
  return x * u;
}

__inline double xldexp(double x, int q) { return ldexpk(x, q); }

__inline int ilogbp1(double d) {
  int m = d < 4.9090934652977266E-91;
  d = m ? 2.037035976334486E90 * d : d;
  int q = (doubleToRawLongBits(d) >> 52) & 0x7ff;
  q = m ? q - (300 + 0x03fe) : q - 0x03fe;
  return q;
}

__inline int xilogb(double d) {
  int e = ilogbp1(xfabs(d)) - 1;
  e = d == 0 ? (-2147483647 - 1) : e;
  e = d == rtengine::RT_INFINITY || d == -rtengine::RT_INFINITY ? 2147483647 : e;
  return e;
}

__inline double upper(double d) {
  return longBitsToDouble(doubleToRawLongBits(d) & 0xfffffffff8000000LL);
}

typedef struct {
  double x, y;
} double2;

typedef struct {
  float x, y;
} float2;

__inline double2 dd(double h, double l) {
  double2 ret;
  ret.x = h; ret.y = l;
  return ret;
}

__inline double2 normalize_d(double2 t) {
  double2 s;

  s.x = t.x + t.y;
  s.y = t.x - s.x + t.y;

  return s;
}

__inline double2 scale_d(double2 d, double s) {
  double2 r;

  r.x = d.x * s;
  r.y = d.y * s;

  return r;
}

__inline double2 add2_ss(double x, double y) {
  double2 r;

  r.x = x + y;
  double v = r.x - x;
  r.y = (x - (r.x - v)) + (y - v);

  return r;
}

__inline double2 add_ds(double2 x, double y) {
  // |x| >= |y|

  double2 r;

  assert(xisnan(x.x) || xisnan(y) || xfabs(x.x) >= xfabs(y));

  r.x = x.x + y;
  r.y = x.x - r.x + y + x.y;

  return r;
}

__inline double2 add2_ds(double2 x, double y) {
  // |x| >= |y|

  double2 r;

  r.x  = x.x + y;
  double v = r.x - x.x;
  r.y = (x.x - (r.x - v)) + (y - v);
  r.y += x.y;

  return r;
}

__inline double2 add_sd(double x, double2 y) {
  // |x| >= |y|

  double2 r;

  assert(xisnan(x) || xisnan(y.x) || xfabs(x) >= xfabs(y.x));

  r.x = x + y.x;
  r.y = x - r.x + y.x + y.y;

  return r;
}

__inline double2 add_dd(double2 x, double2 y) {
  // |x| >= |y|

  double2 r;

  assert(xisnan(x.x) || xisnan(y.x) || xfabs(x.x) >= xfabs(y.x));

  r.x = x.x + y.x;
  r.y = x.x - r.x + y.x + x.y + y.y;

  return r;
}

__inline double2 add2_dd(double2 x, double2 y) {
  double2 r;

  r.x  = x.x + y.x;
  double v = r.x - x.x;
  r.y = (x.x - (r.x - v)) + (y.x - v);
  r.y += x.y + y.y;

  return r;
}

__inline double2 div_dd(double2 n, double2 d) {
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

__inline double2 mul_ss(double x, double y) {
  double xh = upper(x), xl = x - xh;
  double yh = upper(y), yl = y - yh;
  double2 r;

  r.x = x * y;
  r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl;

  return r;
}

__inline double2 mul_ds(double2 x, double y) {
  double xh = upper(x.x), xl = x.x - xh;
  double yh = upper(y  ), yl = y   - yh;
  double2 r;

  r.x = x.x * y;
  r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl + x.y * y;

  return r;
}

__inline double2 mul_dd(double2 x, double2 y) {
  double xh = upper(x.x), xl = x.x - xh;
  double yh = upper(y.x), yl = y.x - yh;
  double2 r;

  r.x = x.x * y.x;
  r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl + x.x * y.y + x.y * y.x;

  return r;
}

__inline double2 squ_d(double2 x) {
  double xh = upper(x.x), xl = x.x - xh;
  double2 r;

  r.x = x.x * x.x;
  r.y = xh * xh - r.x + (xh + xh) * xl + xl * xl + x.x * (x.y + x.y);

  return r;
}

__inline double2 rec_s(double d) {
  double t = 1.0 / d;
  double dh = upper(d), dl = d - dh;
  double th = upper(t), tl = t - th;
  double2 q;

  q.x = t;
  q.y = t * (1 - dh * th - dh * tl - dl * th - dl * tl);

  return q;
}

__inline double2 sqrt_d(double2 d) {
  double t = sqrt(d.x + d.y);
  return scale_d(mul_dd(add2_dd(d, mul_ss(t, t)), rec_s(t)), 0.5);
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
  return mulsign(atan2k(sqrt((1+d)*(1-d)), xfabs(d)), d) + (d < 0 ? rtengine::RT_PI : 0);
}

__inline double xatan(double s) {
  double t, u;
  int q = 0;

  if (s < 0) { s = -s; q = 2; }
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

__inline double xsin(double d) {
  int q;
  double u, s;

  q = (int)xrint(d * rtengine::RT_1_PI);

  d = mla(q, -PI4_A*4, d);
  d = mla(q, -PI4_B*4, d);
  d = mla(q, -PI4_C*4, d);

  s = d * d;

  if ((q & 1) != 0) d = -d;

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

  return u;
}

__inline double xcos(double d) {
  int q;
  double u, s;

  q = 1 + 2*(int)xrint(d * rtengine::RT_1_PI - 0.5);

  d = mla(q, -PI4_A*2, d);
  d = mla(q, -PI4_B*2, d);
  d = mla(q, -PI4_C*2, d);

  s = d * d;

  if ((q & 2) == 0) d = -d;

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

  return u;
}

__inline double2 xsincos(double d) {
  int q;
  double u, s, t;
  double2 r;

  q = (int)xrint(d * (2 * rtengine::RT_1_PI));

  s = d;

  s = mla(-q, PI4_A*2, s);
  s = mla(-q, PI4_B*2, s);
  s = mla(-q, PI4_C*2, s);

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

  u = -1.13615350239097429531523e-11;
  u = mla(u, s, 2.08757471207040055479366e-09);
  u = mla(u, s, -2.75573144028847567498567e-07);
  u = mla(u, s, 2.48015872890001867311915e-05);
  u = mla(u, s, -0.00138888888888714019282329);
  u = mla(u, s, 0.0416666666666665519592062);
  u = mla(u, s, -0.5);

  r.y = u * s + 1;

  if ((q & 1) != 0) { s = r.y; r.y = r.x; r.x = s; }
  if ((q & 2) != 0) { r.x = -r.x; }
  if (((q+1) & 2) != 0) { r.y = -r.y; }

  if (xisinf(d)) { r.x = r.y = rtengine::RT_NAN; }

  return r;
}

__inline double xtan(double d) {
  int q;
  double u, s, x;

  q = (int)xrint(d * (2 * rtengine::RT_1_PI));

  x = mla(q, -PI4_A*2, d);
  x = mla(q, -PI4_B*2, x);
  x = mla(q, -PI4_C*2, x);

  s = x * x;

  if ((q & 1) != 0) x = -x;

  u = 1.01419718511083373224408e-05;
  u = mla(u, s, -2.59519791585924697698614e-05);
  u = mla(u, s, 5.23388081915899855325186e-05);
  u = mla(u, s, -3.05033014433946488225616e-05);
  u = mla(u, s, 7.14707504084242744267497e-05);
  u = mla(u, s, 8.09674518280159187045078e-05);
  u = mla(u, s, 0.000244884931879331847054404);
  u = mla(u, s, 0.000588505168743587154904506);
  u = mla(u, s, 0.00145612788922812427978848);
  u = mla(u, s, 0.00359208743836906619142924);
  u = mla(u, s, 0.00886323944362401618113356);
  u = mla(u, s, 0.0218694882853846389592078);
  u = mla(u, s, 0.0539682539781298417636002);
  u = mla(u, s, 0.133333333333125941821962);
  u = mla(u, s, 0.333333333333334980164153);

  u = mla(s, u * x, x);

  if ((q & 1) != 0) u = 1.0 / u;

  if (xisinf(d)) u = rtengine::RT_NAN;

  return u;
}

__inline double xlog(double d) {
  double x, x2, t, m;
  int e;

  e = ilogbp1(d * 0.7071);
  m = ldexpk(d, -e);

  x = (m-1) / (m+1);
  x2 = x * x;

  t = 0.148197055177935105296783;
  t = mla(t, x2, 0.153108178020442575739679);
  t = mla(t, x2, 0.181837339521549679055568);
  t = mla(t, x2, 0.22222194152736701733275);
  t = mla(t, x2, 0.285714288030134544449368);
  t = mla(t, x2, 0.399999999989941956712869);
  t = mla(t, x2, 0.666666666666685503450651);
  t = mla(t, x2, 2);

  x = x * t + 0.693147180559945286226764 * e;

  if (xisinf(d)) x = rtengine::RT_INFINITY;
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

  if (xisminf(d)) u = 0;

  return u;
}

__inline double2 logk(double d) {
  double2 x, x2;
  double m, t;
  int e;

  e = ilogbp1(d * 0.7071);
  m = ldexpk(d, -e);

  x = div_dd(add2_ss(-1, m), add2_ss(1, m));
  x2 = squ_d(x);

  t = 0.134601987501262130076155;
  t = mla(t, x2.x, 0.132248509032032670243288);
  t = mla(t, x2.x, 0.153883458318096079652524);
  t = mla(t, x2.x, 0.181817427573705403298686);
  t = mla(t, x2.x, 0.222222231326187414840781);
  t = mla(t, x2.x, 0.285714285651261412873718);
  t = mla(t, x2.x, 0.400000000000222439910458);
  t = mla(t, x2.x, 0.666666666666666371239645);

  return add2_dd(mul_ds(dd(0.693147180559945286226764, 2.319046813846299558417771e-17), e),
		 add2_dd(scale_d(x, 2), mul_ds(mul_dd(x2, x), t)));
}

__inline double expk(double2 d) {
  int q = (int)rint((d.x + d.y) * R_LN2);
  double2 s, t;
  double u;

  s = add2_ds(d, q * -L2U);
  s = add2_ds(s, q * -L2L);

  s = normalize_d(s);

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

  t = add_dd(s, mul_ds(squ_d(s), u));

  t = add_sd(1, t);
  return ldexpk(t.x + t.y, q);
}

__inline double xpow(double x, double y) {
  int yisint = (int)y == y;
  int yisodd = (1 & (int)y) != 0 && yisint;

  double result = expk(mul_ds(logk(xfabs(x)), y));

  result = xisnan(result) ? rtengine::RT_INFINITY : result;
  result *=  (x >= 0 ? 1 : (!yisint ? rtengine::RT_NAN : (yisodd ? -1 : 1)));

  double efx = mulsign(xfabs(x) - 1, y);
  if (xisinf(y)) result = efx < 0 ? 0.0 : (efx == 0 ? 1.0 : rtengine::RT_INFINITY);
  if (xisinf(x) || x == 0) result = (yisodd ? sign(x) : 1) * ((x == 0 ? -y : y) < 0 ? 0 : rtengine::RT_INFINITY);
  if (xisnan(x) || xisnan(y)) result = rtengine::RT_NAN;
  if (y == 0 || x == 1) result = 1;

  return result;
}

__inline double2 expk2(double2 d) {
  int q = (int)rint((d.x + d.y) * R_LN2);
  double2 s, t;
  double u;

  s = add2_ds(d, q * -L2U);
  s = add2_ds(s, q * -L2L);

  s = normalize_d(s);

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

  t = add_dd(s, mul_ds(squ_d(s), u));

  t = add_sd(1, t);
  return dd(ldexpk(t.x, q), ldexpk(t.y, q));
}

__inline double xsinh(double x) {
  double y = xfabs(x);
  double2 d = expk2(dd(y, 0));
  d = add2_dd(d, div_dd(dd(-1, 0), d));
  y = (d.x + d.y) * 0.5;

  y = xisinf(x) || xisnan(y) ? rtengine::RT_INFINITY : y;
  y = mulsign(y, x);
  y = xisnan(x) ? rtengine::RT_NAN : y;

  return y;
}

__inline double xcosh(double x) {
  double2 d = expk2(dd(x, 0));
  d = add2_dd(d, div_dd(dd(1, 0), d));
  double y = (d.x + d.y) * 0.5;

  y = xisinf(x) || xisnan(y) ? rtengine::RT_INFINITY : y;
  y = xisnan(x) ? rtengine::RT_NAN : y;

  return y;
}

__inline double xtanh(double x) {
  double y = xfabs(x);
  double2 d = expk2(dd(y, 0));
  double2 e = div_dd(dd(1, 0), d);
  d = div_dd(add2_dd(d, scale_d(e, -1)), add2_dd(d, e));
  y = d.x + d.y;

  y = xisinf(x) || xisnan(y) ? 1.0 : y;
  y = mulsign(y, x);
  y = xisnan(x) ? rtengine::RT_NAN : y;

  return y;
}

__inline double2 logk2(double2 d) {
  double2 x, x2, m;
  double t;
  int e;

  d = normalize_d(d);
  e = ilogbp1(d.x * 0.7071);
  m = scale_d(d, ldexpk(1, -e));

  x = div_dd(add2_ds(m, -1), add2_ds(m, 1));
  x2 = squ_d(x);

  t = 0.134601987501262130076155;
  t = mla(t, x2.x, 0.132248509032032670243288);
  t = mla(t, x2.x, 0.153883458318096079652524);
  t = mla(t, x2.x, 0.181817427573705403298686);
  t = mla(t, x2.x, 0.222222231326187414840781);
  t = mla(t, x2.x, 0.285714285651261412873718);
  t = mla(t, x2.x, 0.400000000000222439910458);
  t = mla(t, x2.x, 0.666666666666666371239645);

  return add2_dd(mul_ds(dd(0.693147180559945286226764, 2.319046813846299558417771e-17), e),
		 add2_dd(scale_d(x, 2), mul_ds(mul_dd(x2, x), t)));
}

__inline double xasinh(double x) {
  double y = xfabs(x);
  double2 d = logk2(add2_ds(sqrt_d(add2_ds(mul_ss(y, y),  1)), y));
  y = d.x + d.y;

  y = xisinf(x) || xisnan(y) ? rtengine::RT_INFINITY : y;
  y = mulsign(y, x);
  y = xisnan(x) ? rtengine::RT_NAN : y;

  return y;
}

__inline double xacosh(double x) {
  double2 d = logk2(add2_ds(sqrt_d(add2_ds(mul_ss(x, x), -1)), x));
  double y = d.x + d.y;

  y = xisinf(x) || xisnan(y) ? rtengine::RT_INFINITY : y;
  y = x == 1.0 ? 0.0 : y;
  y = x < 1.0 ? rtengine::RT_NAN : y;
  y = xisnan(x) ? rtengine::RT_NAN : y;

  return y;
}

__inline double xatanh(double x) {
  double y = xfabs(x);
  double2 d = logk2(div_dd(add2_ss(1, y), add2_ss(1, -y)));
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

  e = ilogbp1(d);
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

__inline double xexp2(double a) {
  double u = expk(mul_ds(dd(0.69314718055994528623, 2.3190468138462995584e-17), a));
  if (xispinf(a)) u = rtengine::RT_INFINITY;
  if (xisminf(a)) u = 0;
  return u;
}

__inline double xexp10(double a) {
  double u = expk(mul_ds(dd(2.3025850929940459011, -2.1707562233822493508e-16), a));
  if (xispinf(a)) u = rtengine::RT_INFINITY;
  if (xisminf(a)) u = 0;
  return u;
}

__inline double xexpm1(double a) {
  double2 d = add2_ds(expk2(dd(a, 0)), -1.0);
  double x = d.x + d.y;
  if (xispinf(a)) x = rtengine::RT_INFINITY;
  if (xisminf(a)) x = -1;
  return x;
}

__inline double xlog10(double a) {
  double2 d = mul_dd(logk(a), dd(0.43429448190325176116, 6.6494347733425473126e-17));
  double x = d.x + d.y;

  if (xisinf(a)) x = rtengine::RT_INFINITY;
  if (a < 0) x = rtengine::RT_NAN;
  if (a == 0) x = -rtengine::RT_INFINITY;

  return x;
}

__inline double xlog1p(double a) {
  double2 d = logk2(add2_ss(a, 1));
  double x = d.x + d.y;

  if (xisinf(a)) x = rtengine::RT_INFINITY;
  if (a < -1) x = rtengine::RT_NAN;
  if (a == -1) x = -rtengine::RT_INFINITY;

  return x;
}

///////////////////////////////////////////

#define PI4_Af 0.78515625f
#define PI4_Bf 0.00024127960205078125f
#define PI4_Cf 6.3329935073852539062e-07f
#define PI4_Df 4.9604681473525147339e-10f

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
__inline float xrintf(float x) { return x < 0 ? (int)(x - 0.5f) : (int)(x + 0.5f); }

__inline int xisnanf(float x) { return x != x; }
__inline int xisinff(float x) { return x == rtengine::RT_INFINITY_F || x == -rtengine::RT_INFINITY_F; }
__inline int xisminff(float x) { return x == -rtengine::RT_INFINITY_F; }
__inline int xispinff(float x) { return x == rtengine::RT_INFINITY_F; }

__inline int ilogbp1f(float d) {
  int m = d < 5.421010862427522E-20f;
  d = m ? 1.8446744073709552E19f * d : d;
  int q = (floatToRawIntBits(d) >> 23) & 0xff;
  q = m ? q - (64 + 0x7e) : q - 0x7e;
  return q;
}

__inline float ldexpkf(float x, int q) {
  float u;
  int m;
  m = q >> 31;
  m = (((m + q) >> 6) - m) << 4;
  q = q - (m << 2);
  u = intBitsToFloat(((int32_t)(m + 0x7f)) << 23);
  u = u * u;
  x = x * u * u;
  u = intBitsToFloat(((int32_t)(q + 0x7f)) << 23);
  return x * u;
}

__inline float xcbrtf(float d) { // max error : 2 ulps
  float x, y, q = 1.0f;
  int e, r;

  e = ilogbp1f(d);
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
  float u, s;

  q = (int)xrintf(d * (float)rtengine::RT_1_PI);

  d = mlaf(q, -PI4_Af*4, d);
  d = mlaf(q, -PI4_Bf*4, d);
  d = mlaf(q, -PI4_Cf*4, d);
  d = mlaf(q, -PI4_Df*4, d);

  s = d * d;

  if ((q & 1) != 0) d = -d;

  u = 2.6083159809786593541503e-06f;
  u = mlaf(u, s, -0.0001981069071916863322258f);
  u = mlaf(u, s, 0.00833307858556509017944336f);
  u = mlaf(u, s, -0.166666597127914428710938f);

  u = mlaf(s, u * d, d);

  return u;
}

__inline float xcosf(float d) {
  int q;
  float u, s;

  q = 1 + 2*(int)xrintf(d * (float)rtengine::RT_1_PI - 0.5f);

  d = mlaf(q, -PI4_Af*2, d);
  d = mlaf(q, -PI4_Bf*2, d);
  d = mlaf(q, -PI4_Cf*2, d);
  d = mlaf(q, -PI4_Df*2, d);

  s = d * d;

  if ((q & 2) == 0) d = -d;

  u = 2.6083159809786593541503e-06f;
  u = mlaf(u, s, -0.0001981069071916863322258f);
  u = mlaf(u, s, 0.00833307858556509017944336f);
  u = mlaf(u, s, -0.166666597127914428710938f);

  u = mlaf(s, u * d, d);

  return u;
}

__inline float2 xsincosf(float d) {
  int q;
  float u, s, t;
  float2 r;

  q = (int)rint(d * ((float)(2 * rtengine::RT_1_PI)));

  s = d;

  s = mlaf(q, -PI4_Af*2, s);
  s = mlaf(q, -PI4_Bf*2, s);
  s = mlaf(q, -PI4_Cf*2, s);
  s = mlaf(q, -PI4_Df*2, s);

  t = s;

  s = s * s;

  u = -0.000195169282960705459117889f;
  u = mlaf(u, s, 0.00833215750753879547119141f);
  u = mlaf(u, s, -0.166666537523269653320312f);
  u = u * s * t;

  r.x = t + u;

  u = -2.71811842367242206819355e-07f;
  u = mlaf(u, s, 2.47990446951007470488548e-05f);
  u = mlaf(u, s, -0.00138888787478208541870117f);
  u = mlaf(u, s, 0.0416666641831398010253906f);
  u = mlaf(u, s, -0.5f);

  r.y = u * s + 1;

  if ((q & 1) != 0) { s = r.y; r.y = r.x; r.x = s; }
  if ((q & 2) != 0) { r.x = -r.x; }
  if (((q+1) & 2) != 0) { r.y = -r.y; }

  if (xisinff(d)) { r.x = r.y = rtengine::RT_NAN_F; }

  return r;
}

__inline float xtanf(float d) {
  int q;
  float u, s, x;

  q = (int)xrintf(d * (float)(2 * rtengine::RT_1_PI));

  x = d;

  x = mlaf(q, -PI4_Af*2, x);
  x = mlaf(q, -PI4_Bf*2, x);
  x = mlaf(q, -PI4_Cf*2, x);
  x = mlaf(q, -PI4_Df*2, x);

  s = x * x;

  if ((q & 1) != 0) x = -x;

  u = 0.00927245803177356719970703f;
  u = mlaf(u, s, 0.00331984995864331722259521f);
  u = mlaf(u, s, 0.0242998078465461730957031f);
  u = mlaf(u, s, 0.0534495301544666290283203f);
  u = mlaf(u, s, 0.133383005857467651367188f);
  u = mlaf(u, s, 0.333331853151321411132812f);

  u = mlaf(s, u * x, x);

  if ((q & 1) != 0) u = 1.0f / u;

  if (xisinff(d)) u = rtengine::RT_NAN_F;

  return u;
}

__inline float xatanf(float s) {
  float t, u;
  int q = 0;

  if (s < 0) { s = -s; q = 2; }
  if (s > 1) { s = 1.0f / s; q |= 1; }

  t = s * s;

  u = 0.00282363896258175373077393f;
  u = mlaf(u, t, -0.0159569028764963150024414f);
  u = mlaf(u, t, 0.0425049886107444763183594f);
  u = mlaf(u, t, -0.0748900920152664184570312f);
  u = mlaf(u, t, 0.106347933411598205566406f);
  u = mlaf(u, t, -0.142027363181114196777344f);
  u = mlaf(u, t, 0.199926957488059997558594f);
  u = mlaf(u, t, -0.333331018686294555664062f);

  t = s + s * (t * u);

  if ((q & 1) != 0) t = 1.570796326794896557998982f - t;
  if ((q & 2) != 0) t = -t;

  return t;
}

__inline float atan2kf(float y, float x) {
  float s, t, u;
  float q = 0.f;

  if (x < 0) { x = -x; q = -2.f; }
  if (y > x) { t = x; x = y; y = -t; q += 1.f; }

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

  t = u * t;
  t = mlaf(t,s,s);
  return mlaf(q,(float)(rtengine::RT_PI_F_2),t);
}

__inline float xatan2f(float y, float x) {
  float r = atan2kf(xfabsf(y), x);

  r = mulsignf(r, x);
  if (xisinff(x) || x == 0) r = rtengine::RT_PI_F/2 - (xisinff(x) ? (signf(x) * (float)(rtengine::RT_PI_F*.5f)) : 0);
  if (xisinff(y)          ) r = rtengine::RT_PI_F/2 - (xisinff(x) ? (signf(x) * (float)(rtengine::RT_PI_F*.25f)) : 0);
  if (              y == 0) r = (signf(x) == -1 ? rtengine::RT_PI_F : 0);

  return xisnanf(x) || xisnanf(y) ? rtengine::RT_NAN_F : mulsignf(r, y);
}

__inline float xasinf(float d) {
  return mulsignf(atan2kf(fabsf(d), sqrtf((1.0f+d)*(1.0f-d))), d);
}

__inline float xacosf(float d) {
  return mulsignf(atan2kf(sqrtf((1.0f+d)*(1.0f-d)), fabsf(d)), d) + (d < 0 ? (float)rtengine::RT_PI : 0.0f);
}

__inline float xlogf(float d) {
  float x, x2, t, m;
  int e;

  e = ilogbp1f(d * 0.7071f);
  m = ldexpkf(d, -e);

  x = (m-1.0f) / (m+1.0f);
  x2 = x * x;

  t = 0.2371599674224853515625f;
  t = mlaf(t, x2, 0.285279005765914916992188f);
  t = mlaf(t, x2, 0.400005519390106201171875f);
  t = mlaf(t, x2, 0.666666567325592041015625f);
  t = mlaf(t, x2, 2.0f);

  x = x * t + 0.693147180559945286226764f * e;

  if (xisinff(d)) x = rtengine::RT_INFINITY_F;
  if (d < 0) x = rtengine::RT_NAN_F;
  if (d == 0) x = -rtengine::RT_INFINITY_F;

  return x;
}

__inline float xexpf(float d) {
  if(d<=-104.0f) return 0.0f;

  int q = (int)xrintf(d * R_LN2f);
  float s, u;

  s = mlaf(q, -L2Uf, d);
  s = mlaf(q, -L2Lf, s);

  u = 0.00136324646882712841033936f;
  u = mlaf(u, s, 0.00836596917361021041870117f);
  u = mlaf(u, s, 0.0416710823774337768554688f);
  u = mlaf(u, s, 0.166665524244308471679688f);
  u = mlaf(u, s, 0.499999850988388061523438f);

  u = mlaf( s, mlaf(s,u,1.f),1.f);
  return ldexpkf(u, q);

}

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



#endif
