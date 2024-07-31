#ifndef __APPROX_H__
#define __APPROX_H__

int firstApprox(unsigned int n, double *X, double *F, double *A, double *Amen);
double firstApprox_func(double x, double a, double b, unsigned int n, double *X, double *A);
int secondApprox(unsigned int n, double *X, double *F, double *A, double *Amen);
double secondApprox_func(double x, double a, double b, unsigned int n, double *X, double *A);

#endif