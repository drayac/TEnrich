#define constant_e (2.71828)
#define PI (3.14159265359)
#define M_2PI  6.283185307179586476925286766559 /* 2*pi */
//#undef min
//#define min(a,b) ((a < b)?a:b)
//#undef max
//#define max(a,b) ((a > b)?a:b)
#define ML_POSINF (1.0 / 0.0)
#define ML_NEGINF ((-1.0) / 0.0)
#define ML_NAN  (0.0 / 0.0)
#ifndef M_LOG10_2
#define M_LOG10_2 0.301029995663981195213738894724 /* log10(2) */
#endif
#define R_FINITE(x)    isfinite(x)
#define R_D__0 (log_p ? ML_NEGINF : 0.)  /* 0 */
#define R_D__1 (log_p ? 0. : 1.)   /* 1 */
#define R_D_exp(x) (log_p ?  (x)  : exp(x)) /* exp(x) */

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E  2.718281828459045235360287471353 /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E  1.442695040888963407359924681002 /* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E 0.434294481903251827651128918917 /* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2  0.693147180559945309417232121458 /* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10  2.302585092994045684017991454684 /* ln(10) */
#endif

#ifndef M_PI
#define M_PI  3.141592653589793238462643383280 /* pi */
#endif

#ifndef M_2PI
#define M_2PI  6.283185307179586476925286766559 /* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2  1.570796326794896619231321691640 /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4  0.785398163397448309615660845820 /* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI  0.318309886183790671537767526745 /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI  0.636619772367581343075535053490 /* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.128379167095512573896158903122 /* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2  1.414213562373095048801688724210 /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362105 /* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3 1.732050807568877293527446341506 /* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32 5.656854249492380195206754896838 /* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2 0.301029995663981195213738894724 /* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI 1.772453850905516027298167483341 /* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI 0.797884560802865355879892119869 /* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI 1.837877066409345483560659472811 /* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI 0.572364942924700087071713675677 /* log(sqrt(pi))
                                                                    == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /* log(sqrt(2*pi))
                                                                   == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2 0.225791352644727432363097614947 /* log(sqrt(pi/2))
                                                                      == log(pi/2)/2 */
#endif

/* FUNCTION HEADERS  */
double logspace_add (double, double) ;

// From Rmath
void bgrat(double, double, double, double, double *, double, int *, bool);
double grat_r(double, double, double, double);
double apser(double, double, double, double);
double bpser(double, double, double, double, bool);
double basym(double, double, double, double, bool);
double fpser(double, double, double, double, bool);
double bup(double, double, double, double, int, double, bool);
double exparg(int);
double psi(double);
double gam1(double);
double gamln1(double);
double betaln(double, double);
double algdiv(double, double);
double brcmp1(int, double, double, double, double, bool);
double brcomp(double, double, double, double, bool);
double rlog1(double);
double bcorr(double, double);
double gamln(double);
double alnrel(double);
double esum(int, double, bool);
double erf__(double);
double rexpm1(double);
double erfc1(int, double);
double gsumln(double, double);

double pbinom(double, double, double, std::string) ;
double pbeta(double, double, double, int, bool) ;
double pbeta_raw(double, double, double, int, bool) ;


