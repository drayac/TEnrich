#include <cmath>
#include <stdlib.h>

#define constant_e (2.71828)
#define PI (3.14159265359)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

#ifndef KFUNC_H_
#define KFUNC_H_
#endif /* KFUNC_H_ */

typedef struct {
    long long n11, n1_, n_1, n;
    double p;
} hgacc_t;

// regularized lower incomplete gamma function, by series expansion
double kf_lgamma(double z);
double kf_erfc(double x);
double _kf_gammap(double s, double z);
double _kf_gammaq(double s, double z);
double kf_gammap(double s, double z);
double kf_gammaq(double s, double z);
double kf_betai_aux(double a, double b, double x);
double kf_betai(double a, double b, double x);
double lbinom(long long n, long long k);
double hypergeo(long long n11, long long n1_, long long n_1, long long n);

// incremental version of hypergeometric distribution
double hypergeo_acc(long long n11, long long n1_, long long n_1, long long n, hgacc_t *aux);
double kt_fisher_exact(long long n11, long long n12, long long n21, long long n22, double *_left, double *_right, double *two);
double log_gamma_f(double) ;
double stirl_fact(double) ;
long double fisher_pval(int a, int b, int c, int d) ;
long double fisher_exact(long long a, long long b, long long c, long long d, std::string) ;

