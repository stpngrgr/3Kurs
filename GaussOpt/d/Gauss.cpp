#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <sys/types.h>
#include <unistd.h>
#include <sys/sysinfo.h>
#include <sched.h>

#include"Gauss.h"





static inline double blockNorm(double* M, int N, int m, double* AM)
{
	double ret = 0.0;
	int k = N / m, l = N - k * m;
	double* block;

	for(int j = 0; j < k; j++)
	{
		for(int jj = 0; jj < m; jj++)
			AM[jj] = 0.;
		block = M + j * m * m;
//		printf("%10.3e -- ", block[0]);
		for(int i = 0; i < k; i++)
		{
			for(int jj = 0; jj < m; jj++)
			{
				for(int ii = 0; ii < m; ii++)
				{
					AM[jj] += fabs(block[ii * m + jj]);
				}
			}
//			printf("%lf ", AM[0]);
			block += N * m;
		}
		if(l != 0)
		{
			block -= j * (m - l) * m;
//			printf("(");
			for(int jj = 0; jj < m; jj++)
			{
				for(int ii = 0; ii < l; ii++)
				{
//					printf("%10.3e ", block[ii * m + jj]);
					AM[jj] += fabs(block[ii * m + jj]);
				}
			}
//			printf(")\n");
//			printf("%lf -- ", AM[0]);
		}
		for(int jj = 0; jj < m; jj++)
		{
			ret = (ret < AM[jj] ? AM[jj] : ret);
		}
//		printf("%d -- %lf\n", j, ret);
//		printf("\n");
	}
	if(l != 0)
	{
		for(int jj = 0; jj < l; jj++)
			AM[jj] = 0.;
		block = M + k * m * m;
		for(int i = 0; i < k; i++)
		{
			for(int jj = 0; jj < l; jj++)
			{
				for(int ii = 0; ii < m; ii++)
				{
					AM[jj] += fabs(block[ii * l + jj]);
				}
			}
			block += N * m;
		}
		block -= k * (m - l) * m;
		for(int jj = 0; jj < l; jj++)
		{
			for(int ii = 0; ii < l; ii++)
			{
				AM[jj] += fabs(block[ii * l + jj]);
			}
		}
		for(int jj = 0; jj < l; jj++)
		{
//			printf("%lf ", AM[jj]);
			ret = (ret < AM[jj] ? AM[jj] : ret);
		}
//		printf("\n");
	}
//	printf("%lf\n", ret);
	return ret;
}

int display(int l, int n, int m, double* A, int r)
{
	int kw = n / m;
	int lw = n - kw * m;
	int kh = l / m;
	int lh = l - kh * m;
	int h = (l > r ? r : l), w = (n > r ? r : n);
	double p;

	printf("\n\n");
	for(int i = 0; i < h; i++)
	{
		for(int j = 0; j < w; j++)
		{
			p = A[(i / m) * n * m + (j / m) * ((i - kh * m) < 0 ? m : lh) * m + (i % m) * ((j - kw * m) < 0 ? m : lw) + (j % m)];
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

static inline double calculateResidual(double* A, double* X, int N, int m)
{
	int k = N / m;
	int l = N - k * m;
	double a = 0.0, b = 0.0, temp = 0.0;

	for(int i = 0; i < N; i++)
	{
		temp = 0.0;
		for(int j = 0; j < N; j++)
			temp += A[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)] * X[j];
		temp -= A[N * N + i];
		a += fabs(temp);
		b += fabs(A[N * N + i]);
	}
	return (a/b);
}

static inline double calculateError(double* X, int N)
{
	double res = 0.0;

	for(int i = 0; i < N; i++)
		res += fabs(X[i] - ((i & 1) == 0 ? 1.0 : 0.0));
	return res/(double)((N+1)/2);
}





int main(int argc, char **argv)
{
	int N, m, r, s, err;
	double* A = nullptr, *X = nullptr, *Am = nullptr, clk, res, mnorm;
	cpu_set_t cpu;
	if((argc != 6) && (argc != 5))
	{
		printf("Bad argc\n Usage: %s N m r s filename\n", argv[0]);
		return -1;
	}
	if((!sscanf(argv[1], "%d", &N)) || (N <= 0))
	{
		printf("Bad N: %s\n Usage: %s N m r s filename\n", argv[1], argv[0]);
		return -1;
	}
	if((!sscanf(argv[2], "%d", &m)) || (m <= 0) || (m > N))
	{
		printf("Bad m: %s\n Usage: %s N m r s filename\n", argv[2], argv[0]);
		return -1;
	}
	else if ((m % 3) != 0)
	{
		printf("m shold be divisible by 3: %s\n Usage: %s N m r s filename\n", argv[2], argv[0]);
		return -1;
	}
	if((!sscanf(argv[3], "%d", &r)) || (r <= 0))
	{
		printf("Bad r: %s\n Usage: %s N m r s filename\n", argv[3], argv[0]);
		return -1;
	}
	if((!sscanf(argv[4], "%d", &s)) || (s < 0))
	{
		printf("Bad s: %s\n Usage: %s N m r s filename\n", argv[4], argv[0]);
		return -1;
	}
	if((s == 0) && (argc != 6))
	{
		printf("Filename missing\n Usage: %s N m r s filename\n", argv[0]);
		return -1;
	}
	CPU_ZERO(&cpu);
	CPU_SET(get_nprocs() - 1, &cpu);
	sched_setaffinity(getpid(), sizeof(cpu), &cpu);
	A = (double*)malloc((N * N + N) * sizeof(double));
	if(nullptr == A)
	{
		printf("A allocation failed\n");
		return -1;
	}
	err = init(A, N, m, s, (s == 0 ? argv[5] : nullptr));
	if(0 != err)
	{
		switch(err)
		{
		case -1:
			printf("A is NULL\n");
			break;
		case -30:
			printf("Can't open file\n");
			break;
		case -31:
			printf("Bad element encountered\n");
			break;
		case -32:
			printf("Not enough elements\n");
			break;
		default:
			printf("Unknown matrix initialization error\n");
		}
		free(A);
		return -1;
	}
	X = (double*)malloc(N * sizeof(double));
	if(nullptr == X)
	{
		printf("X allocation failed\n");
		free(A);
		return -1;
	}
	Am = (double*)malloc(m * m * 3 * sizeof(double));
	if(nullptr == Am)
	{
		printf("Additional memory allocation failed\n");
		free(A);
		free(X);
		return -1;
	}
	display(N, N, m, A, r);

	mnorm = blockNorm(A, N, m, Am);
	clk = clock();
	err = solve(N, m, A, A + (N * N), X, Am, mnorm);
	clk = clock() - clk;
	if(0 != err)
	{
		printf("Algorithm is unapplicable\nFailure on %d step\n", -err);
		free(A);
		free(X);
		free(Am);
		return -1;
	}
	err = init(A, N, m, s, (s == 0 ? argv[5] : nullptr));
	if(0 != err)
	{
		switch(err)
		{
		case -1:
			printf("A is NULL\n");
			break;
		case -30:
			printf("Can't open file\n");
			break;
		case -31:
			printf("Bad element encountered\n");
			break;
		case -32:
			printf("Not enough elements\n");
			break;
		default:
			printf("Unknown matrix initialization error\n");
		}
		free(A);
		free(X);
		free(Am);
		return -1;
	}
	display(N, 1, m, X, r);
	printf("Error : %e\n", calculateError(X, N));
	clk /= CLOCKS_PER_SEC;
	res = calculateResidual(A, X, N, m);
	printf ("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d\n", argv[0], res, clk, s, N, m);
	free(A);
	free(X);
	free(Am);
	return 0;
}


// void mcpy(void* dest, void* from, int size)
// {
// 	char *d = (char*)dest, *f = (char*)from;


// 	for(int i = 0; i < size; i++)
// 	{
// 		d[i] = f[i];
// 	}
// }





//3mult
/*int matrixMult(double* A, double* B, double* C, int p, int q, int r)
{

	for(int i = 0; i < p; i +=3)
	{
		for(int j = 0; j < r; j += 3)
		{
			C[i * q + j] = 0.;
			C[i * q + j + 1] = 0.;
			C[i * q + j + 2] = 0.;
			C[(i + 1) * q + j] = 0.;
			C[(i + 1) * q + j + 1] = 0.;
			C[(i + 1) * q + j + 2] = 0.;
			C[(i + 2) * q + j] = 0.;
			C[(i + 2) * q + j + 1] = 0.;
			C[(i + 2) * q + j + 2] = 0.;
			for(int k = 0; k < q; k++)
			{
				C[i * q + j] += A[i * q + k] * B[k * r + j];
				C[i * q + j + 1] += A[i * q + k] * B[k * r + j + 1];
				C[i * q + j + 2] += A[i * q + k] * B[k * r + j + 2];
				C[(i + 1) * q + j] += A[(i + 1) * q + k] * B[k * r + j];
				C[(i + 1) * q + j + 1] += A[(i + 1) * q + k] * B[k * r + j + 1];
				C[(i + 1) * q + j + 2] += A[(i + 1) * q + k] * B[k * r + j + 2];
				C[(i + 2) * q + j] += A[(i + 2) * q + k] * B[k * r + j];
				C[(i + 2) * q + j + 1] += A[(i + 2) * q + k] * B[k * r + j + 1];
				C[(i + 2) * q + j + 2] += A[(i + 2) * q + k] * B[k * r + j + 2];
			}
		}
	}
	return 0;
}*/


//3mult
/*int matrixMultSub(double* A, double* B, double* C, int p, int q, int r)
{
	double temp[9];

	for(int i = 0; i < p; i += 3)
		for(int j = 0; j < r; j += 3)
		{
			temp[0] = 0.;
			temp[1] = 0.;
			temp[2] = 0.;
			temp[3] = 0.;
			temp[4] = 0.;
			temp[5] = 0.;
			temp[6] = 0.;
			temp[7] = 0.;
			temp[8] = 0.;
			for(int k = 0; k < q; k++)
			{
				temp[0] += A[i * q + k] * B[k * r + j];
				temp[1] += A[i * q + k] * B[k * r + j + 1];
				temp[2] += A[i * q + k] * B[k * r + j + 2];
				temp[3] += A[(i + 1) * q + k] * B[k * r + j];
				temp[4] += A[(i + 1) * q + k] * B[k * r + j + 1];
				temp[5] += A[(i + 1) * q + k] * B[k * r + j + 2];
				temp[6] += A[(i + 2) * q + k] * B[k * r + j];
				temp[7] += A[(i + 2) * q + k] * B[k * r + j + 1];
				temp[8] += A[(i + 2) * q + k] * B[k * r + j + 2];
			}
			C[i * r + j] -= temp[0];
			C[i * r + j + 1] -= temp[1];
			C[i * r + j + 2] -= temp[2];
			C[(i + 1) * r + j] -= temp[3];
			C[(i + 1) * r + j + 1] -= temp[4];
			C[(i + 1) * r + j + 2] -= temp[5];
			C[(i + 2) * r + j] -= temp[6];
			C[(i + 2) * r + j + 1] -= temp[7];
			C[(i + 2) * r + j + 2] -= temp[8];
		}
	return 0;
}*/
