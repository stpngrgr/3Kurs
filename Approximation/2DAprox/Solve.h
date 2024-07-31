#ifndef __SOLVE_H__
#define __SOLVE_H__
#include <pthread.h>
#include <stdio.h>







static inline unsigned int tRowsBeg(unsigned int n, int p, int t)
{
    return (n * t) / p;
}


static inline unsigned int tRowsEnd(unsigned int n, int p, int t)
{
    return (n * (t + 1)) / p;
}


int tSolveRestart(double *A, unsigned int *I, double *b, double *x, double *r, double *u,
                  double *v, double *amt, double *amg, unsigned int n, double eps, int maxit,
                  int p, int t, pthread_barrier_t* barrier);



#endif