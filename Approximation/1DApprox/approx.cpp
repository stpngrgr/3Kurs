#include <stdio.h>
#include <string.h>
#include <math.h>

#include "approx.h"

#define __K__ 4
#define eps 1.e-1

double getSDerivative(int side, unsigned int n, double *X, double *F, double *Amen);
int trid_solve(unsigned int n, double *A, double *B, double *X);

int firstApprox(unsigned int n, double *X, double *F, double *A, double *)
{
    double den = 1.;

    for (unsigned int i = 0; i < n; i++)
    {
        den = 1.;
        for (unsigned int j = 0; j < n; j ++)
        {
            if(i != j)
                den *= (X[i] - X[j]);
        }
        A[i] = F[i] / den;
    }

    return 0;
}

double firstApprox_func(double x, double , double , unsigned int n, double *X, double *A)
{
    double varphi = 1.;
    double sum = 0.;
    unsigned int i;

    for (i = 0; i < n; i++)
    {
        if (fabs(x - X[i]) < 1.e-16)
            break;
        varphi *= (x - X[i]);
        sum += A[i] / (x - X[i]);
    }
    if (i < n)
    {
        // printf("AHA : %10.6e\n", x);
        sum = A[i];
        varphi = 1.;
        for (unsigned int j = 0; j < n; j++)
        {
            if (j != i)
                varphi *= (x - X[j]);
        }
    }
    // printf("(%lf, %lf)\n", x, varphi * sum);

    return varphi * sum;
}


int secondApprox(unsigned int n, double *X, double *F, double *A, double *Amen)
{
    double xi_1, xi_2, xi_3;
    double xi_0 = X[0] - (X[1] - X[0]) / 2, xi_n1 = X[n - 1] + (X[n - 1] - X[n - 2]) / 2;
    double p_0, p_n;
    unsigned int k = (n > __K__ ? __K__ : n);
    // unsigned int k1 = (n > 8 ? 8 : n);
    int err;

    // printf("Enter\n");

    if (k < 2)
        k = 2;

    p_0 = getSDerivative(-1, k, X, F, Amen);
    p_n = getSDerivative(1, k, X + n - k, F + n - k, Amen);
    // p_n = getSDerivative(1, k1, X + n - k1, F + n - k1, Amen);
    // printf("p_0 = %10.16e, p_n = %10.16e\n", p_0, p_n);

    //first additional equation
    xi_1 = (X[1] + X[0]) / 2;
    A[0] = 0.;
    A[1] = 2. / ((xi_1 - xi_0) * (X[0] - xi_0));
    A[2] = 2. / ((xi_1 - xi_0) * (xi_1 - X[0]));
    Amen[0] = p_0 + ((2. * F[0]) / (xi_1 - xi_0)) * (1. / (X[0] - xi_0) + 1. / (xi_1 - X[0]));
    //end
    //
    // xi_1 = (X[1] + X[0]) / 2;
    // A[0] = 0.;
    // A[1] = 2. / ((xi_1 - xi_0) * (X[0] - xi_0));
    // A[2] = 2. / ((xi_1 - xi_0) * (xi_1 - X[0]));
    // Amen[0] = p_0 + ((2. * F[0])/ (xi_1 - xi_0)) * (1. / (X[0] - xi_0) + 1. / (xi_1 - X[0]));
    //


    xi_2 = xi_0;
    xi_3 = (X[1] + X[0]) / 2;
    for (unsigned int i = 1; i < n - 1; i++)
    {
        xi_1 = xi_2;
        xi_2 = xi_3;
        xi_3 = (X[i + 1] + X[i]) / 2;
        A[3 * i] = 1. / (X[i - 1] - xi_1) - 1. / (xi_2 - xi_1);
        A[3 * i + 1] = 1. / (xi_2 - X[i - 1]) + 1. / (xi_2 - xi_1) + 1. / (X[i] - xi_2) + 1. / (xi_3 - xi_2);
        A[3 * i + 2] = 1. / (xi_3 - X[i]) - 1. / (xi_3 - xi_2);
        Amen[i] = (1. / (X[i - 1] - xi_1) + 1. / (xi_2 - X[i - 1])) * F[i - 1] +
            (1. / (X[i] - xi_2) + 1. / (xi_3 - X[i])) * F[i]; 
    }
    xi_1 = xi_2;
    xi_2 = xi_3;
    xi_3 = xi_n1;
    A[3 * (n - 1)] = 1. / (X[(n - 1) - 1] - xi_1) - 1. / (xi_2 - xi_1);
    A[3 * (n - 1) + 1] = 1. / (xi_2 - X[(n - 1) - 1]) + 1. / (xi_2 - xi_1) + 1. / (X[(n - 1)] - xi_2) + 1. / (xi_3 - xi_2);
    A[3 * (n - 1) + 2] = 1. / (xi_3 - X[(n - 1)]) - 1. / (xi_3 - xi_2);
    Amen[(n - 1)] = (1. / (X[(n - 1) - 1] - xi_1) + 1. / (xi_2 - X[(n - 1) - 1])) * F[(n - 1) - 1] +
        (1. / (X[(n - 1)] - xi_2) + 1. / (xi_3 - X[(n - 1)])) * F[(n - 1)]; 


    //second additional equation
    xi_3 = (X[n - 2] + X[n - 1]) / 2;
    A[3 * n] = 2. / ((xi_n1 - xi_3) * (X[n - 1] - xi_3));
    A[3 * n + 1] = 2. / ((xi_n1 - xi_3) * (xi_n1 - X[n - 1]));
    A[3 * n + 2] = 0.;
    Amen[n] = p_n + ((2. * F[n - 1]) / (xi_n1 - xi_3)) * (1. / (X[n - 1] - xi_3) + 1 / (xi_n1 - X[n - 1]));
    //end
    //
    // xi_3 = (X[n - 2] + X[n - 1]) / 2;
    // A[3 * n] = 2. / ((xi_n1 - xi_3) * (X[n - 1] - xi_3));
    // A[3 * n + 1] = 2. / ((xi_n1 - xi_3) * (xi_n1 - X[n - 1]));
    // A[3 * n + 2] = 0.;
    // Amen[n] = p_n + ((2. * F[n - 1])/ (xi_n1 - xi_3)) * (1. / (X[n - 1] - xi_3) + 1 / (xi_n1 - X[n - 1]));
    //


    //debug
    // printf("\n\n e, X:\n");
    // xi_1 = xi_0;
    // for (unsigned int i = 0; i < n - 1; i++)
    // {
    //     printf("%10.16e %10.16e\n", xi_1, X[i]);
    //     xi_1 = (X[i + 1] + X[i]) / 2;
    // }
    // printf("%10.16e %10.16e\n", xi_1, X[n - 1]);
    // xi_1 = xi_n1;
    // printf("%10.16e\n", xi_1);
    // printf("\n\n");
    //


    //debug
    // printf("\n\nA:\n");
    // for (unsigned int i = 0; i < n + 1; i++)
    // {
    //     printf("%10.16e  %10.16e  %10.16e | %10.16e\n", A[3 * i], A[3 * i + 1], A[3 * i + 2], Amen[i]);
    // }
    // printf("\n\n");
    //


    err = trid_solve(n + 1, A, Amen, Amen + n + 1);
    if (err)
        return err;
    
    //debug
    // printf("\n\nX:\n");
    // for (unsigned int i = 0; i < n + 1; i++)
    // {
    //     printf("%10.16e\n", Amen[n + 1 + i]);
    // }
    // printf("\n\n");
    //

    xi_1 = xi_0;
    xi_2 = (X[1] + X[0]) / 2;
    A[0] = Amen[n + 1];
    A[1] = (F[0] - Amen[n + 1]) / (X[0] - xi_1) - ((X[0] - xi_1) / (xi_2 - xi_1)) * 
        ((Amen[n + 2] - F[0]) / (xi_2 - X[0]) - (F[0] - Amen[n + 1]) / (X[0] - xi_1));
    A[2] = (1. / (xi_2 - xi_1)) * ((Amen[n + 2] - F[0]) / (xi_2 - X[0]) - 
        (F[0] - Amen[n + 1]) / (X[0] - xi_1));
    for (unsigned int i = 1; i < n - 1; i++)
    {
        xi_1 = xi_2;
        xi_2 = (X[i + 1] + X[i]) / 2;
        A[3 * i] = Amen[n + 1 + i];
        A[3 * i + 1] = (F[i] - Amen[n + 1 + i]) / (X[i] - xi_1) - ((X[i] - xi_1) / (xi_2 - xi_1)) * 
            ((Amen[n + 2 + i] - F[i]) / (xi_2 - X[i]) - (F[i] - Amen[n + 1 + i]) / (X[i] - xi_1));
        A[3 * i + 2] = (1. / (xi_2 - xi_1)) * ((Amen[n + 2 + i] - F[i]) / (xi_2 - X[i]) - 
            (F[i] - Amen[n + 1 + i]) / (X[i] - xi_1));
    }
    xi_1 = xi_2;
    xi_2 = xi_n1;
    A[3 * (n - 1)] = Amen[2 * n];
    A[3 * (n - 1) + 1] = (F[n - 1] - Amen[2 * n]) / (X[n - 1] - xi_1) - ((X[n - 1] - xi_1) / (xi_2 - xi_1)) * 
        ((Amen[2 * n + 1] - F[n - 1]) / (xi_2 - X[n - 1]) - (F[n - 1] - Amen[2 * n]) / (X[n - 1] - xi_1));
    A[3 * (n - 1) + 2] = (1. / (xi_2 - xi_1)) * ((Amen[2 * n + 1] - F[n - 1]) / (xi_2 - X[n - 1]) - 
        (F[n - 1] - Amen[2 * n]) / (X[n - 1] - xi_1));

    //debug
    // printf("\n\nCof:\n");
    // for (unsigned int i = 0; i < n; i++)
    // {
    //     printf("%10.16e  %10.16e  %10.16e\n", A[3 * i], A[3 * i + 1], A[3 * i + 2]);
    // }
    // printf("\n\n");
    //

    // printf("Leave\n");

    return 0;
}


double secondApprox_func(double x, double , double , unsigned int n, double *X, double *A)
{
    unsigned int i, left = 0, right = n - 1;
    double xi, ret;

    while (right > left + 1)
    {
        i = (right + left) / 2;
        if (x < X[i])
            right = i;
        else
            left = i;
    }
    i = right;

    xi = (X[i - 1] + X[i]) / 2;

    if (x < xi)
    {
        if (i > 1)
            xi = (X[i - 2] + X[i - 1]) / 2;
        else
            xi = X[0] - (X[1] - X[0]) / 2;
        ret = A[3 * (i - 1)] + A[3 * (i - 1) + 1] * (x - xi) + A[3 * (i - 1) + 2] * (x - xi) * (x - xi);  
    }
    else
    {
        ret = A[3 * i] + A[3 * i + 1] * (x - xi) + A[3 * i + 2] * (x - xi) * (x - xi);
    }
    
    return ret;
}





int trid_solve(unsigned int n, double *A, double *B, double *X)
{
    double norm = fabs(A[1]) + fabs(A[3]), temp;

    //debug
    // for (int i = 0; i < n; i++)
    // {
    //     // printf("%lf  %lf  %lf\n", A[3 * i], A[3 * i + 1], A[3 * i + 2]);
    //     printf("%lf\n", B[i]);
    // }

    for (unsigned int i = 1; i < n - 1; i++)
    {
        temp = fabs(A[3 * (i - 1) + 2]) + fabs(A[3 * i + 1]) + fabs(A[3 * (i + 1)]);
        norm = (temp > norm ? temp : norm);
    }
    temp = fabs(A[3 * (n - 2) + 2]) + fabs(A[3 * (n - 1) + 1]);
    norm = (temp > norm ? temp : norm);

    memcpy(X, B, n * sizeof(double));
    for (unsigned int i = 0; i < n - 1; i++)
    {
        if (fabs(A[3 * i + 1]) < 1.e-16 * norm)
            return -i - 1;
        A[3 * i + 2] /= A[3 * i + 1];
        X[i] /= A[3 * i + 1];
        A[3 * (i + 1) + 1] -= A[3 * i + 2] * A[3 * (i + 1)];
        X[i + 1] -= X[i] * A[3 * (i + 1)];
    }
    if (fabs(A[3 * (n - 1) + 1]) < 1.e-16 * norm)
        return -n;
    X[n - 1] /= A[3 * (n - 1) + 1];
    for (int i = n - 1; i > 0; i--)
    {
        X[i - 1] -= X[i] * A[3 * (i - 1) + 2];
    }

    return 0;
}


double getDerivative(double x, unsigned int n, double *X, double *, double *Amen)
{ 
    double y_0 = 0., y_1 = 0., y_2 = 0.;
    double y0 = 0., y1 = 0., y2 = 0.;
    double m_0 = 0., m_1 = 0., m_2 = 0.;
    double ul1 = 0., ul2 = 0.;
    double ret = 0.;

    y_0 = firstApprox_func(x - eps, ul1, ul2, n, X, Amen);
    y_1 = firstApprox_func(x - eps * 2., ul1, ul2, n, X, Amen);
    y_2 = firstApprox_func(x - eps * 3., ul1, ul2, n, X, Amen);
    y0 = firstApprox_func(x + eps, ul1, ul2, n, X, Amen);
    y1 = firstApprox_func(x + eps * 2., ul1, ul2, n, X, Amen);
    y2 = firstApprox_func(x + eps * 3., ul1, ul2, n, X, Amen);

    m_0 = (y0 - y_0) / 2.;
    m_1 = (y1 - y_1) / 4.;
    m_2 = (y2 - y_2) / 6.;

    ret = ((m_0 * 15.) - (m_1 * 6.) + m_2) / (10. * eps);

    //debug
    // printf("\nFirst %.5e : %10.16e %10.16e %10.16e %10.16e %10.16e %10.16e \n %10.16e\n\n", x, y_2, y_1, y_0, y0, y1, y2, ret);
    //

    return ret;
}


double getSDerivative(int side, unsigned int n, double *X, double *F, double *Amen)
{
    double y_0 = 0., y_1 = 0., y_2 = 0.;
    double y0 = 0., y1 = 0., y2 = 0.;
    double m_0 = 0., m_1 = 0., m_2 = 0.;
    double ret = 0.;

    firstApprox(n, X, F, Amen, nullptr);
    if (side <= 0)
    {
        y0 = getDerivative(X[0] + eps, n, X, F, Amen);
        y1 = getDerivative(X[0] + eps * 2., n, X, F, Amen);
        y2 = getDerivative(X[0] + eps * 3., n, X, F, Amen);
        y_0 = getDerivative(X[0] - eps, n, X, F, Amen);
        y_1 = getDerivative(X[0] - eps * 2., n, X, F, Amen);
        y_2 = getDerivative(X[0] - eps * 3., n, X, F, Amen);

        m_0 = (y0 - y_0) / 2.;
        m_1 = (y1 - y_1) / 4.;
        m_2 = (y2 - y_2) / 6.;

        ret = ((m_0 * 15.) - (m_1 * 6.) + m_2) / (10. * eps);

        //debug
        // printf("\nSecond left %10.16e\n\n", ret);
        // printf("%10.16e %10.16e %10.16e %10.16e %10.16e %10.16e \n\n", y_2, y_1, y_0, y0, y1, y2);
        //
    }
    else
    {
        y0 = getDerivative(X[n - 1] + eps, n, X, F, Amen);
        y1 = getDerivative(X[n - 1] + eps * 2., n, X, F, Amen);
        y2 = getDerivative(X[n - 1] + eps * 3., n, X, F, Amen);
        y_0 = getDerivative(X[n - 1] - eps, n, X, F, Amen);
        y_1 = getDerivative(X[n - 1] - eps * 2., n, X, F, Amen);
        y_2 = getDerivative(X[n - 1] - eps * 3., n, X, F, Amen);

        m_0 = (y0 - y_0) / 2.;
        m_1 = (y1 - y_1) / 4.;
        m_2 = (y2 - y_2) / 6.;

        ret = ((m_0 * 15.) - (m_1 * 6.) + m_2) / (10. * eps);

        //debug
        // printf("\nSecond right %10.16e\n\n", ret);
        // printf("%10.16e %10.16e %10.16e %10.16e %10.16e %10.16e \n\n", y_2, y_1, y_0, y0, y1, y2);
        //
    }

    return ret;
}