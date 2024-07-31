#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include <sys/types.h>
#include <unistd.h>
#include <sys/sysinfo.h>
#include <sched.h>

#ifdef __DEBUG__
#include <stddef.h>
#include <fenv.h>
#include <unistd.h>
#endif

#include"Eigenvalues.h"

//int display(int l, int N, int m, double* A, int r);
double calculateResidual(double* A, double* X, int N, int m);
double calculateError(double* X, int N);




int display(int l, int n, double* A, int r)
{
	int h = (l > r ? r : l);
	int w = (n > r ? r : n);
	double p;

	printf("\n\n");
	for(int i = 0; i < h; i++)
	{
		for(int j = 0; j < w; j++)
		{
			p = A[i * n + j];
			//p = (fabs(p) > 1e-16 ? p : 0.0);
			if(0 >= printf(" %10.3e", p))
				return -1;
		}
		if(0 >= printf("\n"))
			return -1;
	}
	printf("\n\n");
	return 0;
}




int calculateResidual(double* arr, int len, double* res1, double* res2, double cof)
{
	double ret1 = *res1, ret2 = *res2, temp = 0.;
	for(int i = 0; i < len; i++)
	{
		ret1 -= arr[i];
		temp += arr[i] * arr[i];
	}
	ret2 -= sqrt(temp);
	*res1 = fabs(ret1) / cof;
	*res2 = fabs(ret2) / cof;
	return 0;
}




int main(int argc, char **argv)
{
	int N, r, s, err;
	double *A = nullptr, *Am = nullptr, e;
	double res1 = 0., res2 = 0., cof = 0.;
	double tridclk, solveclk;
	int iter;
	node* list;
//	cpu_set_t cpu;
	
	#ifdef __DEBUG__
		feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
	#endif
	
	if((argc != 6) && (argc != 5))
	{
		printf("ERROR : Bad argc\n Usage: %s N m e k filename\n", argv[0]);
		return -1;
	}
	if((!sscanf(argv[1], "%d", &N)) || (N <= 0))
	{
		printf("ERROR : Bad N: %s\n Usage: %s N m e k filename\n", argv[1], argv[0]);
		return -1;
	}
	if((!sscanf(argv[2], "%d", &r)) || (r <= 0))
	{
		printf("ERROR : Bad m: %s\n Usage: %s N m e k filename\n", argv[2], argv[0]);
		return -1;
	}
	if((!sscanf(argv[3], "%lf", &e)) || (e <= 0) || (fabs(e) < e_mash))
	{
		printf("ERROR : Bad e: %s\n Usage: %s N m e k filename\n", argv[3], argv[0]);
		return -1;
	}
	if((!sscanf(argv[4], "%d", &s)) || (s < 0))
	{
		printf("ERROR : Bad k: %s\n Usage: %s N m e k filename\n", argv[4], argv[0]);
		return -1;
	}
	if((s == 0) && (argc != 6))
	{
		printf("ERROR : Filename missing\n Usage: %s N m e k filename\n", argv[0]);
		return -1;
	}
//	CPU_ZERO(&cpu);
//	CPU_SET(get_nprocs() - 1, &cpu);
//	sched_setaffinity(getpid(), sizeof(cpu), &cpu);
	A = new double[N * N];
	if(nullptr == A)
	{
		printf("ERROR : A allocation failed\n");
		return -1;
	}
	err = init(A, N, &res1, &res2, &cof, s, (s == 0 ? argv[5] : nullptr));
	if(0 != err)
	{
		switch(err)
		{
		case -1:
			printf("ERROR : A is NULL\n");
			break;
		case -2:
			printf("ERROR : Algorithm is unapplicable\nMatrix is non symmetric\n");
			break;
		case -30:
			printf("ERROR : Can't open file\n");
			break;
		case -31:
			printf("ERROR : Bad element encountered\n");
			break;
		case -32:
			printf("ERROR : Not enough elements\n");
			break;
		default:
			printf("ERROR : Unknown matrix initialization error\n");
		}
		delete[] A;
		return -1;
	}
	Am = new double[3 * N];
	if(nullptr == Am)
	{
		printf("ERROR : Am allocation failed\n");
		delete[] A;
		return -1;
	}
	list = new node[N];
	if(nullptr == list)
	{
		printf("ERROR : Am allocation failed\n");
		delete[] A;
		delete[] Am;
		return -1;
	}
	display(N, N, A, r);
	
	tridclk = clock();
	err = tridiagonalize(A, N, Am, Am + N, Am + 2 * N, cof);
	tridclk = clock() - tridclk;
	if (0 != err)
	{
		printf("ERROR : What?\n");
		delete[] A;
		delete[] Am;
		delete[] list;
		return -1;
	}
	solveclk = clock();
	iter = SolveRec(A, N, e, Am, list, cof);
	solveclk = clock() - solveclk;
	if(-1 == iter)
	{
		printf("ERROR : Whats taking it so long?\n");
		delete[] A;
		delete[] Am;
		delete[] list;
		return -1;
	}
	if(-2 == iter)
	{
		printf("ERROR : Pardonte?\n");
		delete[] A;
		delete[] Am;
		delete[] list;
		return -1;
	}
//	err = init(A, N, m, s, (s == 0 ? argv[5] : nullptr));
//	if(0 != err)
//	{
//		switch(err)
//		{
//		case -1:
//			printf("A is NULL\n");
//			break;
//		case -30:
//			printf("Can't open file\n");
//			break;
//		case -31:
//			printf("Bad element encountered\n");
//			break;
//		case -32:
//			printf("Not enough elements\n");
//			break;
//		default:
//			printf("Unknown matrix initialization error\n");
//		}
//		free(A);
//		free(Am);
//		return -1;
//	}
//	printf("Error : %e\n", calculateError(X, N));
	display(1, N, Am, r);
	tridclk /= CLOCKS_PER_SEC;
	solveclk /= CLOCKS_PER_SEC;
	calculateResidual(Am, N, &res1, &res2, cof);
	printf("%s: Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], res1, res2, iter, iter / N, tridclk, solveclk);
//	printf ("%s: residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", argv[0], res, clk, s, N, m);
	delete[] A;
	delete[] Am;
	delete[] list;
	return 0;
}





//double calculateResidual(double* A, double* X, int N, int m)
//{
//	int k = N / m;
//	int l = N - k * m;
//	double a = 0.0, b = 0.0, temp = 0.0;
//
//	for(int i = 0; i < N; i++)
//	{
//		temp = 0.0;
//		for(int j = 0; j < N; j++)
//			temp += A[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)] * X[j];
//		temp -= A[N * N + i];
//		a += fabs(temp);
//		b += fabs(A[N * N + i]);
//	}
//	return (a/b);
//}
//
//double calculateError(double* X, int N)
//{
//	double res = 0.0;
//
//	for(int i = 0; i < N; i++)
//		res += fabs(X[i] - ((i & 1) == 0 ? 1.0 : 0.0));
//	return res/(double)((N+1)/2);
//}
//
//
//double blockNorm(double* M, int N, int m)
//{
//	double ret = 0.0, temp = 0.0;
//	int k = N / m, l = N - k * m;
//	for(int j = 0; j < N; j++)
//	{
//		temp = 0.0;
//		for(int i = 0; i < N; i++)
//		{
//			temp += fabs(M[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)]);
//		}
//		ret = (ret < temp ? temp : ret);
//	}
//	return ret;
//}
