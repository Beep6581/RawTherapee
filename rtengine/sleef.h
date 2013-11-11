typedef struct {
  double x, y;
} double2;

typedef struct {
  float x, y;
} float2;

double xsin(double d);
double xcos(double d);
double2 xsincos(double d);
double xtan(double d);
double xasin(double s);
double xacos(double s);
double xatan(double s);
double xatan2(double y, double x);
double xlog(double d);
double xexp(double d);
double xpow(double x, double y);

double xsinh(double x);
double xcosh(double x);
double xtanh(double x);
double xasinh(double x);
double xacosh(double x);
double xatanh(double x);
double xldexp(double x, int q);
int xilogb(double d);

double xfma(double x, double y, double z);
double xsqrt(double d);
double xcbrt(double d);

double xexp2(double a);
double xexp10(double a);
double xexpm1(double a);
double xlog10(double a);
double xlog1p(double a);

float xsinf(float d);
float xcosf(float d);
float2 xsincosf(float d);
float xtanf(float d);
float xasinf(float s);
float xacosf(float s);
float xatanf(float s);
float xatan2f(float y, float x);
float xlogf(float d);
float xexpf(float d);
float xpowf(float x, float y);
float xcbrtf(float d);
