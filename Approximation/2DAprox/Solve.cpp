#include "Solve.h"







// static void display(double *x, unsigned int n)
// {
//     printf("\n\n ---------------- \n\n");
//     for (unsigned int i = 0; i < n; i++)
//     {
//         printf("%e\n", x[i]);
//     }
//     printf("\n\n ---------------- \n\n");
// }




static inline double reduce_sum_ordered(double *sarr, int p, pthread_barrier_t *barrier)
{
	// static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    double ret;
	int barret;
    int j;

	//lock
	barret = pthread_barrier_wait(barrier);
	if (barret == PTHREAD_BARRIER_SERIAL_THREAD)
    {
        for (int len = p; len > 1; len = (len + 1) / 2)
        {
            for (j = 0; j < len / 2; j++)
            {
                ret = sarr[2 * j] + sarr[2 * j + 1];
                sarr[j] = ret;
            }
            if ((len % 2) != 0)
                sarr[j] = sarr[len - 1];
        }
    }
	pthread_barrier_wait(barrier);
	ret = sarr[0];
	pthread_barrier_wait(barrier);
	return ret;
}





void tMatVecMult(double *A, unsigned int *I, double *x, double *y,
                 unsigned int n, int p, int t, pthread_barrier_t* barrier)
{
    double s;

    for (unsigned int i = tRowsBeg(n, p, t); i < tRowsEnd(n, p, t); i++)
    {
        s = A[i] * x[i];
        for (unsigned int j = I[i]; j < I[i + 1]; j++)
        {
            s += A[j] * x[I[j]];
        }
        y[i] = s;
    }  
    //sync
    pthread_barrier_wait(barrier);
}


void tPrecond(double *A, unsigned int */*I*/, double *x, double *y, 
              unsigned int n, int p, int t, pthread_barrier_t* barrier)
{
    for (unsigned int i = tRowsBeg(n, p, t); i < tRowsEnd(n, p, t); i++)
    {
        x[i] = y[i] / A[i];
    }
    //sync
    pthread_barrier_wait(barrier);
}


void tAdjVec(double *x, double *y, double alpha, unsigned int n, int p, 
    int t, pthread_barrier_t* barrier)
{
    for (unsigned int i = tRowsBeg(n, p, t); i < tRowsEnd(n, p, t); i++)
    {
        x[i] -= alpha * y[i];
    }
    //sync
    pthread_barrier_wait(barrier);
}


double tScalarProd(double *x, double *y, unsigned int n, double *amt, double *amg,
    int p, int t, pthread_barrier_t* barrier)
{
    double ret;
    unsigned int amb = tRowsBeg(n, p, t), i, j = amb;


    // printf("\n\n ----------------------------------- \n\n");
    for (j = 0; j < (tRowsEnd(n, p, t) - tRowsBeg(n, p, t)) / 2; j++)
    {
        // printf("%e -- %e\n%e -- %e\n", x[amb + 2 * j], y[amb + 2 * j], x[amb + 2 * j + 1], y[amb + 2 * j + 1]);
        amt[j] = x[amb + 2 * j] * y[amb + 2 * j] + x[amb + 2 * j + 1] * y[amb + 2 * j + 1];
    }
    if (((tRowsEnd(n, p, t) - tRowsBeg(n, p, t)) % 2) != 0)
        amt[j] = x[tRowsEnd(n, p, t) - 1] * y[tRowsEnd(n, p, t) - 1];

    
    for (i = (tRowsEnd(n, p, t) - tRowsBeg(n, p, t) + 1) / 2; i > 1; i = (i + 1) / 2)
    {
        for (j = 0; j < i / 2; j++)
        {
            ret = amt[2 * j] + amt[2 * j + 1];
            amt[j] = ret;
        }
        if (i % 2 != 0)
            amt[j] = amt[i - 1];
    }
    amg[t] = amt[0];
    //reduce
    ret = reduce_sum_ordered(amg, p, barrier);
    // printf("\n\n ----------------------------------- \n\n");

    return ret;
}




// double tScalarProd(double *x, double *y, unsigned int n, double *amt, double *amg,
//     int p, int t, pthread_barrier_t* barrier)
// {
//     double ret;
//     // unsigned int amb = tRowsBeg(n, p, t), i, j = amb;


//     // printf("\n\n ----------------------------------- \n\n");
//     // for (j = 0; j < (tRowsEnd(n, p, t) - tRowsBeg(n, p, t)) / 2; j++)
//     // {
//     //     // printf("%e -- %e\n%e -- %e\n", x[amb + 2 * j], y[amb + 2 * j], x[amb + 2 * j + 1], y[amb + 2 * j + 1]);
//     //     amt[j] = x[amb + 2 * j] * y[amb + 2 * j] + x[amb + 2 * j + 1] * y[amb + 2 * j + 1];
//     // }
//     // if ((tRowsEnd(n, p, t) % 2) != 0)
//     //     amt[j] = x[tRowsEnd(n, p, t) - 1] * y[tRowsEnd(n, p, t) - 1];

    
//     // for (i = (tRowsEnd(n, p, t) - tRowsBeg(n, p, t) + 1) / 2; i > 1; i = (i + 1) / 2)
//     // {
//     //     for (j = 0; j < i / 2; j++)
//     //     {
//     //         ret = amt[2 * j] + amt[2 * j + 1];
//     //         amt[j] = ret;
//     //     }
//     //     if (i % 2 != 0)
//     //         amt[j] = amt[i - 1];
//     // }

//     amt[0] = 0.;
//     for (unsigned int i = tRowsBeg(n, p, t); i < tRowsEnd(n, p, t); i++)
//     {
//         amt[0] += x[i] * y[i];
//     }
//     amg[t] = amt[0];
//     //reduce
//     ret = reduce_sum_ordered(amg, p, barrier);
//     // printf("\n\n ----------------------------------- \n\n");

//     return ret;
// }





int tSolve(double *A, unsigned int *I, double *b, double *x, double *r, double *u,
           double *v, double *amt, double *amg, unsigned int n, double eps, int maxit,
           int p, int t, pthread_barrier_t* barrier)
{
    double prec = eps * eps * tScalarProd(b, b, n, amt, amg, p, t, barrier);
    double c1 = 0., c2 = 0.;
    double tau = 0.;
    int it;

    // printf("AAAAAAAAAAAA  %e, %e\n", eps, tScalarProd(b, b, n, amt, amg, p, t, barrier));
    tMatVecMult(A, I, x, r, n, p, t, barrier);
    tAdjVec(r, b, 1., n, p, t, barrier);

    for (it = 0; it < maxit; it++)
    {
        // display(r, n);
        // display(A, n);
        tPrecond(A, I, v, r, n, p, t, barrier);
        // display(v, n);
        tMatVecMult(A, I, v, u, n, p, t, barrier);
        // display(u, n);
        c1 = tScalarProd(u, r, n, amt, amg, p, t, barrier);
        c2 = tScalarProd(u, u, n, amt, amg, p, t, barrier);

        // printf("%d %d :::::::::::     %e, %e, %e\n", t, it, c1, c2, prec);
        if ((c1 < prec) && (c2 < prec))
            break;

        tau = c1 / c2;
        tAdjVec(x, v, tau, n, p, t, barrier);
        tAdjVec(r, u, tau, n, p, t, barrier);
    }

    if (it >= maxit)
        return -1;

    return it;
}


int tSolveRestart(double *A, unsigned int *I, double *b, double *x, double *r, double *u,
                  double *v, double *amt, double *amg, unsigned int n, double eps, int maxit,
                  int p, int t, pthread_barrier_t* barrier)
{
    int clock = 10, it;
    int ret;

    for (it = 0; it < maxit; it += clock)
    {
        ret = tSolve(A, I, b, x, r, u, v, amt, amg, n, eps, clock, p, t, barrier);
        if (ret >= 0)
        {
            it += ret;
            break;
        }
    }

    // printf("kkkkkkkkkkkkkkkkkit %d\n", it);

    if (it >= maxit)
        return -1;
    
    return it;
}
