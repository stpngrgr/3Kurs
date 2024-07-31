#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <sys/sysinfo.h>
//#include <sched.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <pthread.h>

#ifdef __DEBUG__
#include <stddef.h>
#include <fenv.h>
#include <unistd.h>
#endif

#include"Gauss.h"

//int display(int l, int N, int m, double* A, int r);
//double calculateError(double* X, int N, int m);

static inline int args_init(thread_args* args, int p, int N, int m, double* A, double* B, double* X, double *F, int s, int r, char* filename, pthread_barrier_t* barrier)
{
	for(int i = 0; i < p; i++)
	{
		args[i].N = N;
		args[i].m = m;
		args[i].A = A;
		args[i].B = B;
		args[i].X = X;
		args[i].F = F;
		args[i].p = p;
		args[i].t = i;
		args[i].s = s;
		args[i].r = r;
		args[i].res = -1.;
		args[i].filename = filename;
		args[i].barrier = barrier;
	}
	return 0;
}

static inline double calculateError(double* X, int N, int m)
{
	int k = N / m, l = N - k * m;
	double res = 0.0;
	double* block;

	for(int i = 0; i < k; i++)
	{
		block = getBlock(0, i, 1, N, m, X);
		for(int j = 0; j < m; j++)
		{
			res += fabs(block[j] - (((i * m + j) & 1) == 0 ? 1.0 : 0.0));
		}
	}
	if(l != 0)
	{
		block = getBlock(0, k, 1, N, m, X);
		for(int j = 0; j < l; j++)
		{
			res += fabs(block[j] - (((k * m + j) & 1) == 0 ? 1.0 : 0.0));
		}
	}
	return res/(double)((N+1)/2);
}

int main(int argc, char **argv)
{
	int N, m, r, s, p, k, l, err;
	double *A = nullptr, *B = nullptr, *X = nullptr, *F = nullptr;
	char* filename = nullptr;

	thread_args* args;
	pthread_t* ids;
	pthread_barrier_t barrier;

#ifdef __DEBUG__
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
#endif

	if((argc != 7) && (argc != 6))
	{
		printf("Bad argc\n Usage: %s N m p r s filename\n", argv[0]);
		return -1;
	}
	if((!sscanf(argv[1], "%d", &N)) || (N <= 0))
	{
		printf("Bad N: %s\n Usage: %s N m p r s filename\n", argv[1], argv[0]);
		return -1;
	}
	if((!sscanf(argv[2], "%d", &m)) || (m <= 0) || (m > N))
	{
		printf("Bad m: %s\n Usage: %s N m p r s filename\n", argv[2], argv[0]);
		return -1;
	}
	else if ((m % 3) != 0)
	{
		printf("m shold be divisible by 3: %s\n Usage: %s N m r s filename\n", argv[2], argv[0]);
		return -1;
	}
	if((!sscanf(argv[3], "%d", &p)) || (p <= 0))
	{
		printf("Bad N: %s\n Usage: %s N m p r s filename\n", argv[3], argv[0]);
		return -1;
	}
	if((!sscanf(argv[4], "%d", &r)) || (r <= 0))
	{
		printf("Bad r: %s\n Usage: %s N m p r s filename\n", argv[4], argv[0]);
		return -1;
	}
	if((!sscanf(argv[5], "%d", &s)) || (s < 0))
	{
		printf("Bad s: %s\n Usage: %s N m p r s filename\n", argv[5], argv[0]);
		return -1;
	}
	if(s == 0)
	{
		if(argc != 7)
		{
			printf("Filename missing\n Usage: %s N m r s filename\n", argv[0]);
			return -1;
		}
		else
		{
			filename = argv[6];
		}
	}
	k = N / m;
	l = N - k * m;
	A = new double[(blockSizeD(m, m) * k + blockSizeD(m, l)) * k + blockSizeD(l, m) * k + blockSizeD(l, l)];
	if(nullptr == A)
	{
		printf("A allocation failed\n");
		return -1;
	}
	B = new double[blockSizeD(1, m) * k + blockSizeD(1, l)];
	if(nullptr == B)
	{
		printf("B allocation failed\n");
		delete[] A;
		return -1;
	}
	X = new double[blockSizeD(1, m) * k + blockSizeD(1, l)];
	if(nullptr == X)
	{
		printf("X allocation failed\n");
		delete[] A;
		delete[] B;
		return -1;
	}
	F = new double[(blockSizeD(m, m) * (k + 2) + blockSizeD(m, l)) * p];
	if(nullptr == F)
	{
		printf("Additional memory allocation failed\n");
		delete[] A;
		delete[] B;
		delete[] X;
		return -1;
	}
	args = new thread_args[p];
	if(nullptr == args)
	{
		printf("failed to allocate memory for args\n");
		delete[] A;
		delete[] B;
		delete[] X;
		delete[] F;
		return -2;
	}
	ids = new pthread_t[p];
	if(nullptr == ids)
	{
		printf("failed to allocate memory for ids\n");
		delete[] A;
		delete[] B;
		delete[] X;
		delete[] F;
		delete[] args;
		return -1;
	}
	pthread_barrier_init(&barrier, 0, p);
	args_init(args, p, N, m, A, B, X, F, s, r, filename, &barrier);
	for(int i = 1; i < p; i++)
	{
		err = pthread_create(ids + i, 0, thread_func, args + i);
		if(err)
		{
			printf("pthread_create failed on step %d: %s\n", i, strerror(err));
			delete[] A;
			delete[] B;
			delete[] X;
			delete[] F;
			delete[] ids;
			delete[] args;
			return -1;
		}
	}
	thread_func(args + 0);
	if(0 != args[0].ret)
	{
		if(args[0].ret < 10)
		{
			switch(args[0].ret)
			{
			case 1:
				printf("A is NULL\n");
				break;
			case 2:
				printf("Can't open file\n");
				break;
			case 3:
				printf("Bad element encountered\n");
				break;
			case 4:
				printf("Not enough elements\n");
				break;
			default:
				printf("Unknown matrix initialization error\n");
			}
		}
		else
		{
			printf("Algorithm is unapplicable\nFailure on %d step\n", args[0].ret - 9);
		}

		delete[] A;
		delete[] B;
		delete[] X;
		delete[] F;
		delete[] ids;
		delete[] args;
		pthread_barrier_destroy(&barrier);
		return -1;
	}

	
	#ifdef __RESO__
	if (p != 1)
		args[0].res = 0.;
	#endif
	
	if (p == 1)
		args[0].res = -1.;

	for(int i = 0; i < p; i++)
	{
		printf("Thread %d time: %.2f\n", i, 0./*args[i].thread_time*/);
	}
	printf("\n\n");
	printf("Error : %e\n", calculateError(X, N, m));
	//ppppppp
	printf ("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d p = %d\n", argv[0], args[0].res, 0./*args[0].full_time*/, s, N, m, p);
	//sleep(3);
	pthread_barrier_destroy(&barrier);
	delete[] A;
	delete[] B;
	delete[] X;
	delete[] F;
	delete[] ids;
	delete[] args;
	return 0;
}




// int display(int l, int n, int m, double* A, int r)
// {
// 	int kw = n / m;
// 	int lw = n - kw * m;
// 	int kh = l / m;
// 	int lh = l - kh * m;
// 	int h = (l > r ? r : l), w = (n > r ? r : n);
// 	double p;

// 	printf("\n\n");
// 	for(int i = 0; i < h; i++)
// 	{
// 		for(int j = 0; j < w; j++)
// 		{
// 			p = A[(i / m) * n * m + (j / m) * ((i - kh * m) < 0 ? m : lh) * m + (i % m) * ((j - kw * m) < 0 ? m : lw) + (j % m)];
// 			//p = (fabs(p) > 1e-16 ? p : 0.0);
// 			if(0 >= printf(" %10.3e", p))
// 				return -1;
// 		}
// 		if(0 >= printf("\n"))
// 			return -1;
// 	}
// 	printf("\n\n");
// 	return 0;
// }

// double calculateResidual(double* A, double* B, double* X, double* F, int N, int m)
// {
// 	int k = N / m;
// 	int l = N - k * m;
// 	double a = 0.0, b = 0.0, temp = 0.0;

// 	for(int i = 0; i < N; i++)
// 	{
// 		temp = 0.0;
// 		for(int j = 0; j < N; j++)
// 			temp += A[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)] * X[j];
// 		temp -= B[i];
// 		a += fabs(temp);
// 		b += fabs(A[N * N + i]);
// 	}
// 	return (a/b);
// }

// double calculateError(double* X, int N, int m)
// {
// 	double res = 0.0;

// 	for(int i = 0; i < N; i++)
// 		res += fabs(X[i] - ((i & 1) == 0 ? 1.0 : 0.0));
// 	return res/(double)((N+1)/2);
// }



int setE(double* M, int m)
{
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < m; j++)
		{
			M[i * m + j] = (i != j ? 0 : 1);
		}
	}
	return 0;
}

int setO(double* M, int n, int m)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++)
		{
			M[i * m + j] = 0;
		}
	}
	return 0;
}





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
