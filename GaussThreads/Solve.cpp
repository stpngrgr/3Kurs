#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "Gauss.h"
#include <string.h>
#include <stdio.h>
#include <sys/sysinfo.h>
#include <sched.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>







static inline int display(int l, int n, int m, double* A, int r)
{
	int kw = n / m;
	int lw = n - kw * m;
	int h = (l > r ? r : l), w = (n > r ? r : n);
	int relh, bh;
	double p;
	double* block;

	printf("\n\n");
	for(int i = 0; i < h; i++)
	{
		bh = (i / m);
		relh = i % m;
		for(int j = 0; j < w; j++)
		{
			block = getBlock(bh, (j / m), l, n, m, A);
			p = block[relh * ((j / m) < kw ? m : lw) + (j % m)];
			if(0 >= printf(" %10.3e", p))
				return -1;

		}
		if(0 >= printf("\n"))
			return -1;
	}
	printf("\n\n");
	return 0;
}

static inline void calculateResidual(double* A, double* B, double* X, double* F, int N, int m, int p, int t, double* a, double* b)
{
	int k = N / m;
	int l = N - k * m;
	double* block, *Xblock, *Fblock;

	*a = 0.;
	*b = 0.;
	Fblock = getBlock(t, 0, m * p, N + 2 * m, m, F);

			//if(t == 0)
				//printf("OOO %lf %lf\n", *a, *b);
	for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
	{
		for(int j = 0; j < m; j++)
		{
			Fblock[j] = 0.;
		}
		for(int j = 0; j < k; j++)
		{
			block = getBlock(t + q * p, j, N, N, m, A);
			Xblock = getBlock(0, j, 1, N, m, X);
			for(int ii = 0; ii < m; ii++)
			{
				for(int jj = 0; jj < m; jj++)
				{
					Fblock[ii] += block[ii * m + jj] * Xblock[jj];
				}
			}
			//if(t == 0)
				//printf("AAA %lf %lf\n", *a, *b);
		}
		if(l != 0)
		{
			block = getBlock(t + q * p, k, N, N, m, A);
			Xblock = getBlock(0, k, 1, N, m, X);
			for(int ii = 0; ii < m; ii++)
			{
				for(int jj = 0; jj < l; jj++)
				{
					Fblock[ii] += block[ii * l + jj] * Xblock[jj];
				}
			}
			//if(t == 0)
				//printf("AAA1 %lf %lf\n", *a, *b);
		}
		Xblock = getBlock(0, t + q * p, 1, N, m, B);
		for(int j = 0; j < m; j++)
		{
			Fblock[j] -= Xblock[j];
			*a += fabs(Fblock[j]);
			*b += fabs(Xblock[j]);
		}
			//if(t == 0)
				//printf("AAA2 %lf %lf\n", *a, *b);
	}
	if((t == (p - 1)) && (l != 0))
	{
		for(int j = 0; j < l; j++)
		{
			Fblock[j] = 0.;
		}
		for(int j = 0; j < k; j++)
		{
			block = getBlock(k, j, N, N, m, A);
			Xblock = getBlock(0, j, 1, N, m, X);
			for(int ii = 0; ii < l; ii++)
			{
				for(int jj = 0; jj < m; jj++)
				{
					Fblock[ii] += block[ii * m + jj] * Xblock[jj];
				}
			}
		}
				//printf("BBB %lf %lf\n", *a, *b);
		block = getBlock(k, k, N, N, m, A);
		Xblock = getBlock(0, k, 1, N, m, X);
		for(int ii = 0; ii < l; ii++)
		{
			for(int jj = 0; jj < l; jj++)
			{
				Fblock[ii] += block[ii * l + jj] * Xblock[jj];
			}
		}
				//printf("BBB1 %lf %lf\n", *a, *b);
		Xblock = getBlock(0, k, 1, N, m, B);
		for(int j = 0; j < l; j++)
		{
			Fblock[j] -= Xblock[j];
			*a += fabs(Fblock[j]);
			*b += fabs(Xblock[j]);
		}
				//printf("BBB2 %lf %lf\n", *a, *b);
	}

	return;
}





static inline int matrixMultRewriteRight(double* A, double* B, double* C, int p, int r)
{
	int prem = p % 3, rrem = r % 3;

	for(int j = 0; j < r - rrem; j += 3)
	{
		for(int i = 0; i < p - prem; i += 3)
		{
			C[i] = 0.;
			C[i + 1] = 0.;
			C[i + 2] = 0.;
			C[p + i] = 0.;
			C[p + i + 1] = 0.;
			C[p + i + 2] = 0.;
			C[2 * p + i] = 0.;
			C[2 * p + i + 1] = 0.;
			C[2 * p + i + 2] = 0.;
			for(int k = 0; k < p; k++)
			{
				C[i] += A[i * p + k] * B[k * r + j];
				C[p + i] += A[i * p + k] * B[k * r + j + 1];
				C[2 * p + i] += A[i * p + k] * B[k * r + j + 2];
				C[i + 1] += A[(i + 1) * p + k] * B[k * r + j];
				C[p + i + 1] += A[(i + 1) * p + k] * B[k * r + j + 1];
				C[2 * p + i + 1] += A[(i + 1) * p + k] * B[k * r + j + 2];
				C[i + 2] += A[(i + 2) * p + k] * B[k * r + j];
				C[p + i + 2] += A[(i + 2) * p + k] * B[k * r + j + 1];
				C[2 * p + i + 2] += A[(i + 2) * p + k] * B[k * r + j + 2];
			}
		}
		for(int i = p - prem; i < p; i++)
		{
			C[i] = 0.;
			C[p + i] = 0.;
			C[2 * p + i] = 0.;
			for(int k = 0; k < p; k++)
			{
				C[i] += A[i * p + k] * B[k * r + j];
				C[p + i] += A[(i + 1) * p + k] * B[k * r + j];
				C[2 * p + i] += A[(i + 2) * p + k] * B[k * r + j];
			}
		}
		for(int i = 0; i < p - prem; i += 3)
		{
			B[i * r + j] = C[i];
			B[i * r + j + 1] = C[p + i];
			B[i * r + j + 2] = C[2 * p + i];
			B[(i + 1) * r + j] = C[i + 1];
			B[(i + 1) * r + j + 1] = C[p + i + 1];
			B[(i + 1) * r + j + 2] = C[2 * p + i + 1];
			B[(i + 2) * r + j] = C[i + 2];
			B[(i + 2) * r + j + 1] = C[p + i + 2];
			B[(i + 2) * r + j + 2] = C[2 * p + i + 2];
		}
		for(int i = p - prem; i < p; i++)
		{
			B[i * r + j] = C[i];
			B[(i + 1) * r + j] = C[p + i];
			B[(i + 2) * r + j] = C[2 * p + i];
		}
	}
	for(int j = r - rrem; j < r; j++)
	{
		for(int i = 0; i < p - prem; i += 3)
		{
			C[i] = 0.;
			C[i + 1] = 0.;
			C[i + 2] = 0.;
			for(int k = 0; k < p; k++)
			{
				C[i] += A[i * p + k] * B[k * r + j];
				C[i + 1] += A[(i + 1) * p + k] * B[k * r + j];
				C[i + 2]+= A[(i + 2) * p + k] * B[k * r + j];
			}
		}
		for(int i = p - prem; i < p; i++)
		{
			C[i] = 0.;
			for(int k = 0; k < p; k++)
			{
				C[i] += A[i * p + k] * B[k * r + j];
			}
		}
		for(int i = 0; i < p - prem; i += 3)
		{
			B[i * r + j] = C[i];
			B[(i + 1) * r + j] = C[i + 1];
			B[(i + 2) * r + j] = C[i + 2];
		}
		for(int i = p - prem; i < p; i++)
		{
			B[i * r + j] = C[i];
		}
	}
	return 0;
}

//nnon3mult
static inline int matrixMultSubFrag(double* A, double* B, double* C, int p, int q, int r)
{
	double temp[9];
	int prem = p % 3, rrem = r % 3;

	for(int i = 0; i < p - prem; i += 3)
	{
		for(int j = 0; j < r - rrem; j += 3)
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
		for(int j = r - rrem; j < r; j++)
		{
			temp[0] = 0.;
			temp[1] = 0.;
			temp[2] = 0.;
			for(int k = 0; k < q; k++)
			{
				temp[0] += A[i * q + k] * B[k * r + j];
				temp[1] += A[(i + 1) * q + k] * B[k * r + j];
				temp[2] += A[(i + 2) * q + k] * B[k * r + j];
			}
			C[i * r + j] -= temp[0];
			C[(i + 1) * r + j] -= temp[1];
			C[(i + 2) * r + j] -= temp[2];
		}
	}
	for(int i = p - prem; i < p; i++)
	{
		for(int j = 0; j < r - rrem; j += 3)
		{
			temp[0] = 0.;
			temp[1] = 0.;
			temp[2] = 0.;
			for(int k = 0; k < q; k++)
			{
				temp[0] += A[i * q + k] * B[k * r + j];
				temp[1] += A[i * q + k] * B[k * r + j + 1];
				temp[2] += A[i * q + k] * B[k * r + j + 2];
			}
			C[i * r + j] -= temp[0];
			C[i * r + j + 1] -= temp[1];
			C[i * r + j + 2] -= temp[2];
		}
		for(int j = r - rrem; j < r; j++)
		{
			temp[0] = 0.;
			for(int k = 0; k < q; k++)
			{
				temp[0] += A[i * q + k] * B[k * r + j];
			}
			C[i * r + j] -= temp[0];
		}
	}
	return 0;
}

static inline int matrixInverse(double* A, double* Ai, int p, double scale)
{
	int mei = -1;
	double maxe = 0.0, temp = 0.0;
	
	for(int i = 0; i < p; i++)
		for(int j = 0; j < p; j++)
			Ai[i * p + j] = (i != j ? 0.0 : 1.0);
	//Down
	for(int i = 0; i < p; i++)
	{
		//display(p, p, p, A, p);
		//display(p, p, p, Ai, p);
		//Mselect
		mei = -1;
		maxe = 0.0;
		for(int k = i; k < p; k++)
		{
			temp = A[k * p + i];
			if(fabs(temp) < (1e-16) * scale)
				continue;
			else if ((mei < 0) || (maxe < fabs(temp)))
			{
				mei = k;
				maxe = temp;
			}	
		}
		if(mei >= 0)
		{
			//RowSwap
			if(mei != i)
			{
				for(int j = 0; j < p; j++)
				{
					temp = A[mei * p + j];
					maxe = Ai[mei * p + j];
					A[mei * p + j] = A[i * p + j];
					Ai[mei * p + j] = Ai[i * p + j];
					A[i * p + j] = temp;
					Ai[i * p + j] = maxe;
				}
			}
			//RowMult
			for(int j = i + 1; j < p; j++)
			{
				A[i * p + j] /= A[i * p + i];
			}
			for(int j = 0; j < p; j++)
				Ai[i * p + j] /= A[i * p + i];
//			A[i * p + i] = 1.0;
			//Substr
			for(int j = i + 1; j < p; j++)
			{
				for(int k = i + 1; k < p; k++)
				{
					A[j * p + k] -= A[j * p + i] * A[i * p + k];
				}
				for(int k = 0; k < p; k++)
					Ai[j * p + k] -= A[j * p + i] * Ai[i * p + k];
//				A[j * p + i] = 0.0;
			}
		}
		else
		{
			return -i - 1;
		}
	}
	
	//Up
	for(int i = p - 1; i > 0; i--)
		for(int j = 0; j < i; j++)
			for(int k = 0; k < p; k++)
				Ai[j * p + k] -= A[j * p + i] * Ai[i * p + k];
	//display(p, p, p, Ai, 10);
	return 0;
}

static inline double blockNorm(double* M, int n, int m)
{
	double ret = 0.0, temp = 0.0;
	for(int j = 0; j < m; j++)
	{
		temp = 0.0;
		for(int i = 0; i < n; i++)
		{
			temp += fabs(M[i * m + j]);
		}
		ret = (ret < temp ? temp : ret);
	}
	return ret;
}

static inline double norm(double* M, int N, int m, double* AM)
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







static inline double thread_time()
{
	struct rusage r;
	getrusage(RUSAGE_THREAD, &r);
	return r.ru_utime.tv_sec + r.ru_utime.tv_usec/(1e+6);
}

static inline double full_time()
{
	struct timeval r;
	gettimeofday(&r, 0);
	return r.tv_sec + r.tv_usec/(1e+6);
}

// static inline int blockSizeB(int rows, int cols)
// {
// 	return blockSizeD(rows, cols) * sizeof(double);
// }

static inline void dcpy(double* dest, double* from, int size)
{
	for(int i = 0; i < size; i++)
	{
		dest[i] = from[i];
	}
}

static inline void mcpy(void* dest, void* from, int size)
{
	char *d = (char*)dest, *f = (char*)from;


	for(int i = 0; i < size; i++)
	{
		d[i] = f[i];
	}
}

static inline void reduce_err(int* ret, pthread_barrier_t* barrier)
{
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static int global = 0;
	int barret;


	//lock
	pthread_mutex_lock(&m);
	if((global == 0) && (*ret != 0))
		global = *ret;
	pthread_mutex_unlock(&m);
	pthread_barrier_wait(barrier);
	*ret = global;
	barret = pthread_barrier_wait(barrier);
	if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		global = 0;
	}
	pthread_barrier_wait(barrier);
	return;
}

static inline void reduce_max(double* s, pthread_barrier_t* barrier)
{
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static double max = 0.;
	static int set = 0;
	int barret;

	//lock
	pthread_mutex_lock(&m);
	if((set == 0) || (max < *s))
	{
		max = *s;
		set = 1;
	}
	pthread_mutex_unlock(&m);
	pthread_barrier_wait(barrier);
	*s = max;
	barret = pthread_barrier_wait(barrier);
	if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		max = 0.;
		set = 0;
	}
	pthread_barrier_wait(barrier);
	return;
}

static inline void reduce_sum_err(double* d, pthread_barrier_t* barrier)
{
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static double sum = 0.;
	int barret;

	//lock
	pthread_mutex_lock(&m);
	if(sum >= 0.)
	{
		if(*d >= 0.)
			sum += *d;
		else
			sum = *d;
	}
	pthread_mutex_unlock(&m);
	pthread_barrier_wait(barrier);
	*d = sum;
	barret = pthread_barrier_wait(barrier);
	if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		sum = 0.;
	}
	pthread_barrier_wait(barrier);
	return;
}

static inline void reduce_min_di(double* d, int* i, pthread_barrier_t* barrier)
{
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static double mind = 0.;
	static int mini = -1;
	int barret;

	pthread_mutex_lock(&m);
	if( (*i >= 0) && ( (mini < 0) || (*d < mind) ) )
	{
		mind = *d;
		mini = *i;
	}
	pthread_mutex_unlock(&m);
	pthread_barrier_wait(barrier);
	*i = mini;
	*d = mind;
	barret = pthread_barrier_wait(barrier);
	if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		mind = 0.;
		mini = -1;
	}
	pthread_barrier_wait(barrier);
	return;
}


static inline void reduce_res(double* a, double* b, pthread_barrier_t* barrier)
{
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static double suma = 0.;
	static double sumb = 0.;
	int barret;

	//lock
	pthread_mutex_lock(&m);
	suma += *a;
	sumb += *b;
	pthread_mutex_unlock(&m);
	barret = pthread_barrier_wait(barrier);
	if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
		suma = suma / sumb;
	pthread_barrier_wait(barrier);
	*a = suma;
	barret = pthread_barrier_wait(barrier);
	if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		suma = 0.;
		sumb = 0.;
	}
	pthread_barrier_wait(barrier);
	return;
}








static inline int findCardinal_threads(int p, int t, double* A, double* B, int N, int m, double mtrNorm, int s, double* F, pthread_barrier_t* barrier)
{
	int k = N / m, l = N - k * m;
	int err = 0, mei = -1;
	double *Ablock, *Ftemp, *Fcopy;
	double min = 0., tempmin = 0.;


	//local
	Ftemp = getBlock(t, k, m * p, N + 2 * m, m, F);
	Fcopy = getBlock(t, k + 1, m * p, N + 2 * m, m, F);

	for(int q = (s / p) + ((s % p) > t ? 1 : 0); ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
	{
		Ablock = getBlock(t + q * p, s, N, N, m, A);
		memcpy(Fcopy, Ablock, m * m * sizeof(double));		
		//debugDelete
		{
			printf("AHA %d\n", q);
			display(m, m, m, Ablock, m);
			display(m, m, m, getBlock(t, 0, m * p, N + 2 * m, m, F), m);
			display(m, m, m, Fcopy, m);
		}
		err = matrixInverse(Fcopy, Ftemp, m, mtrNorm);
		//debugDelete
			printf ("%d\n", err);
			display(m, m, m, Ftemp, m);
		if(err == 0)
		{
			tempmin = blockNorm(Ftemp, m, m);
			if((mei < 0) || (min > tempmin))
			{
				min = tempmin;
				mei = t + q * p;
				Ablock = getBlock(t, 0, m * p, N + 2 * m, m, F);
				memcpy(Ablock, Ftemp, m * m * sizeof(double));
			}
		}
	}
	//sync
	reduce_min_di(&min, &mei, barrier);
	if(mei >= 0)
	{
		if(t == (s % p))
		{
			if(s != mei)
			{
				if(((k - s) - ((k - s) / 2)) > 0)
				{
					Ablock = getBlock(s, s, N, N, m, A);
					Ftemp = getBlock(mei, s, N, N, m, A);
					Fcopy = getBlock(t, 1, m * p, N + 2 * m, m, F);
					memcpy(Fcopy, Ftemp, (blockSizeD(m, m) * ((k - s) - ((k - s) / 2))) * sizeof(double));
					memcpy(Ftemp, Ablock, (blockSizeD(m, m) * ((k - s) - ((k - s) / 2))) * sizeof(double));
					memcpy(Ablock, Fcopy, (blockSizeD(m, m) * ((k - s) - ((k - s) / 2))) * sizeof(double));
				}
				Ablock = getBlock(0, s, 1, N, m, B);
				Ftemp = getBlock(0, mei, 1, N, m, B);
				Fcopy = getBlock(t, 1, m * p, N + 2 * m, m, F);
				memcpy(Fcopy, Ftemp, blockSizeD(m, 1) * sizeof(double));
				memcpy(Ftemp, Ablock, blockSizeD(m, 1) * sizeof(double));
				memcpy(Ablock, Fcopy, blockSizeD(m, 1) * sizeof(double));
			}
			Ftemp = getBlock(t, k, m * p, N + 2 * m, m, F);
			Fcopy = getBlock((mei % p), 0, m * p, N + 2 * m, m, F);
			memcpy(Ftemp, Fcopy, m * m * sizeof(double));
		}
		if(t == (mei % p))
		{
			if((s != mei) && ((l != 0) || (((k - s) / 2) > 0)))
			{
				Ablock = getBlock(s, s + ((k - s) - ((k - s) / 2)), N, N, m, A);
				Ftemp = getBlock(mei, s + ((k - s) - ((k - s) / 2)), N, N, m, A);
				Fcopy = getBlock(t, 1, m * p, N + 2 * m, m, F);
				memcpy(Fcopy, Ftemp, (blockSizeD(m, m) * ((k - s) / 2) + blockSizeD(m, l)) * sizeof(double));
				memcpy(Ftemp, Ablock, (blockSizeD(m, m) * ((k - s) / 2) + blockSizeD(m, l)) * sizeof(double));
				memcpy(Ablock, Fcopy, (blockSizeD(m, m) * ((k - s) / 2) + blockSizeD(m, l)) * sizeof(double));
			}
			if(t != (s % p))
			{
				Ftemp = getBlock(t, k, m * p, N + 2 * m, m, F);
				Fcopy = getBlock((mei % p), 0, m * p, N + 2 * m, m, F);
				memcpy(Ftemp, Fcopy, m * m * sizeof(double));
			}
		}
	}
	pthread_barrier_wait(barrier);
	if(mei >= 0)
	{
		if((t != (s % p)) && (t != (mei % p)))
		{
			Ftemp = getBlock(t, k, m * p, N + 2 * m, m, F);
			Fcopy = getBlock((s % p), k, m * p, N + 2 * m, m, F);
			memcpy(Ftemp, Fcopy, m * m * sizeof(double));
		}
		if(t != (s % p))
		{
			for(int q = (s / p) + ((s % p) + 1 > t ? 1 : 0); ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
			{
				Ftemp = getBlock(t, t + q * p, m * p, N + 2 * m, m, F);
				Fcopy = getBlock(s, t + q * p, N, N, m, A);
				memcpy(Ftemp, Fcopy, m * m * sizeof(double));
			}
			if((t == (p - 1)) && (l != 0))
			{
				Ftemp = getBlock(t, k + 2, m * p, N + 2 * m, m, F);
				Fcopy = getBlock(s, k, N, N, m, A);
				memcpy(Ftemp, Fcopy, m * l * sizeof(double));
			}
		}
		err = 0;
	}
	else
	{
		err = -1;
	}
	
	//debugDelete
	printf("Card : %d\n", mei);
	
	return err;
}

static inline int multRow_threads(int p, int t, double* A, double* B, int N, int m, int s, double* F, pthread_barrier_t* barrier)
{
	int k = N / m, l = N - k * m;
	double *block, *inverse, *temp;


	inverse = getBlock(t, k, m * p, N + 2 * m, m, F);

	//debugDelete
	printf("inv :\n");
	display(m, m, m, inverse, m);


	temp = getBlock(t, k + 1, m * p, N + 2 * m, m, F);
	if(t != (s % p))
	{
		for(int q = (s / p) + ((s % p) + 1 > t ? 1 : 0); ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
		{
			block = getBlock(t, t + q * p, m * p, N + 2 * m, m, F);
			matrixMultRewriteRight(inverse, block, temp, m, m);
		}
		if((t == (p - 1)) && (l != 0))
		{
			block = getBlock(t, k + 2, m * p, N + 2 * m, m, F);
			matrixMultRewriteRight(inverse, block, temp, m, l);
		}
	}
	else
	{
		for(int q = (s / p) + ((s % p) + 1 > t ? 1 : 0); ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
		{
			block = getBlock(s, t + q * p, N, N, m, A);
			matrixMultRewriteRight(inverse, block, temp, m, m);
		}
		if((t == (p - 1)) && (l != 0))
		{
			block = getBlock(s, k, N, N, m, A);
			matrixMultRewriteRight(inverse, block, temp, m, l);
		}
		block = getBlock(0, s, 1, N, m, B);
		matrixMultRewriteRight(inverse, block, temp, m, 1);
	}

	//sync
	pthread_barrier_wait(barrier);

	if(t != (s % p))
	{
		for(int ti = 0; ti < p; ti++)
		{
			if((ti != t) && (ti != (s % p)))
			{
				for(int q = (s / p) + ((s % p) + 1 > ti ? 1 : 0); ((k - ti - 1) >= 0) && (q <= (k - ti - 1) / p); q++)
				{
					block = getBlock(t, ti + q * p, m * p, N + 2 * m, m, F);
					temp = getBlock(ti, ti + q * p, m * p, N + 2 * m, m, F);
					memcpy(block, temp, m * m * sizeof(double));
				}
				if((ti == (p - 1)) && (l != 0))
				{
					block = getBlock(t, k + 2, m * p, N + 2 * m, m, F);
					temp = getBlock(p - 1, k + 2, m * p, N + 2 * m, m, F);
					memcpy(block, temp, m * l * sizeof(double));
				}
			}
			else if(ti == (s % p))
			{
				for(int q = (s / p) + ((s % p) + 1 > ti ? 1 : 0); ((k - ti - 1) >= 0) && (q <= (k - ti - 1) / p); q++)
				{
					block = getBlock(t, ti + q * p, m * p, N + 2 * m, m, F);
					temp = getBlock(s, ti + q * p, N, N, m, A);
					memcpy(block, temp, m * m * sizeof(double));
				}
				if((ti == (p - 1)) && (l != 0))
				{
					block = getBlock(t, k + 2, m * p, N + 2 * m, m, F);
					temp = getBlock(s, k, N, N, m, A);
					memcpy(block, temp, m * l * sizeof(double));
				}
				block = getBlock(t, k, m * p, N + 2 * m, m, F);
				temp = getBlock(0, s, 1, N, m, B);
				memcpy(block, temp, m * sizeof(double));
			}
		}
	}
	else
	{
		for(int ti = 0; ti < p; ti++)
		{
			if(ti != t)
			{
				for(int q = (s / p) + ((s % p) + 1 > ti ? 1 : 0); ((k - ti - 1) >= 0) && (q <= (k - ti - 1) / p); q++)
				{
					block = getBlock(s, ti + q * p, N, N, m, A);
					temp = getBlock(ti, ti + q * p, m * p, N + 2 * m, m, F);
					memcpy(block, temp, m * m * sizeof(double));
				}
				if((ti == (p - 1)) && (l != 0))
				{
					block = getBlock(s, k, N, N, m, A);
					temp = getBlock(p - 1, k + 2, m * p, N + 2 * m, m, F);
					memcpy(block, temp, m * l * sizeof(double));
				}
			}
		}
	}

	return 0;
}

static inline int downStep_threads(double* A, double* B, double* F, int N, int m, int p, int t, int s)
{
	int k = N / m, l = N - k * m;
	double *AblockMult, *AblockSub, *Fblock;


	if(t != (s % p))
	{
		for(int q = (s / p) + (((s % p) + 1) > t ? 1 : 0); ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
		{
			AblockMult = getBlock(t + q * p, s, N, N, m, A);
			for(int j = s + 1; j < k; j++)
			{
				Fblock = getBlock(t, j, m * p, N + 2 * m, m, F);
				AblockSub = getBlock(t + q * p, j, N, N, m, A);
				matrixMultSubFrag(AblockMult, Fblock, AblockSub, m, m, m);
			}
			if(l != 0)
			{
				Fblock = getBlock(t, k + 2, m * p, N + 2 * m, m, F);
				AblockSub = getBlock(t + q * p, k, N, N, m, A);
				matrixMultSubFrag(AblockMult, Fblock, AblockSub, m, m, l);
			}
			Fblock = getBlock(t, k, m * p, N + 2 * m, m, F);
			AblockSub = getBlock(0, t + q * p, 1, N, m, B);
			matrixMultSubFrag(AblockMult, Fblock, AblockSub, m, m, 1);
		}
		if((t == (p - 1)) && (l != 0))
		{
			AblockMult = getBlock(k, s, N, N, m, A);
			for(int j = s + 1; j < k; j++)
			{
				Fblock = getBlock(t, j, m * p, N + 2 * m, m, F);
				AblockSub = getBlock(k, j, N, N, m, A);
				matrixMultSubFrag(AblockMult, Fblock, AblockSub, l, m, m);
			}
			if(l != 0)
			{
				Fblock = getBlock(t, k + 2, m * p, N + 2 * m, m, F);
				AblockSub = getBlock(k, k, N, N, m, A);
				
				#ifdef __DEBUG__
				printf("Step %d :: \n", s);
				display(l, l, l, AblockSub, l);
				#endif
				matrixMultSubFrag(AblockMult, Fblock, AblockSub, l, m, l);
				#ifdef __DEBUG__
				printf("Step %d over :: \n", s);
				display(l, l, l, AblockSub, l);
				#endif
			}
			Fblock = getBlock(t, k, m * p, N + 2 * m, m, F);
			AblockSub = getBlock(0, k, 1, N, m, B);
			matrixMultSubFrag(AblockMult, Fblock, AblockSub, l, m, 1);
		}
	}
	else
	{
		for(int q = (s / p) + (((s % p) + 1) > t ? 1 : 0); ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
		{
			AblockMult = getBlock(t + q * p, s, N, N, m, A);
			for(int j = s + 1; j < k; j++)
			{
				Fblock = getBlock(s, j, N, N, m, A);
				AblockSub = getBlock(t + q * p, j, N, N, m, A);
				matrixMultSubFrag(AblockMult, Fblock, AblockSub, m, m, m);
			}
			if(l != 0)
			{
				Fblock = getBlock(s, k, N, N, m, A);
				AblockSub = getBlock(t + q * p, k, N, N, m, A);
				matrixMultSubFrag(AblockMult, Fblock, AblockSub, m, m, l);
			}
			Fblock = getBlock(0, s, 1, N, m, B);
			AblockSub = getBlock(0, t + q * p, 1, N, m, B);
			matrixMultSubFrag(AblockMult, Fblock, AblockSub, m, m, 1);
		}
		if((t == (p - 1)) && (l != 0))
		{
			AblockMult = getBlock(k, s, N, N, m, A);
			for(int j = s + 1; j < k; j++)
			{
				Fblock = getBlock(s, j, N, N, m, A);
				AblockSub = getBlock(k, j, N, N, m, A);
				matrixMultSubFrag(AblockMult, Fblock, AblockSub, l, m, m);
			}
			if(l != 0)
			{
				Fblock = getBlock(s, k, N, N, m, A);
				AblockSub = getBlock(k, k, N, N, m, A);
				matrixMultSubFrag(AblockMult, Fblock, AblockSub, l, m, l);
			}
			Fblock = getBlock(0, s, 1, N, m, B);
			AblockSub = getBlock(0, k, 1, N, m, B);
			matrixMultSubFrag(AblockMult, Fblock, AblockSub, l, m, 1);
		}
	}
	
	
	//debugDelete	
	printf("B :\n");
	display(N, 1, m, B, N);


	return 0;
}

static inline int lastDownFirstUp_threads(double* A, double* B, int N, int m, int p, int t, double norm, double* F, pthread_barrier_t* barrier)
{
	int k = N / m, l = N - k * m;
	double *block, *inverse, *temp;
	int ret = 0;

	if(t == (p - 1))
	{
		if(l != 0)
		{
			block = getBlock(k, k, N, N, m, A);
			inverse = getBlock(t, k, m * p, N + 2 * m, m, F);
			temp = getBlock(t, k + 1, m * p, N + 2 * m, m, F);
			
			
			#ifdef __DEBUG__
			printf("Last :\n");
			display(l, l, l, block, l);
			#endif

			ret = matrixInverse(block, inverse, l, norm);
			
			
			#ifdef __DEBUG__
			printf("LastInv :\n");
			display(l, l, l, inverse,l);
			#endif
			
			if(ret == 0)
			{
				block = getBlock(0, k, 1, N, m, B);
				matrixMultRewriteRight(inverse, block, temp, l, 1);
			}
		}
	}
	reduce_err(&ret, barrier);
	if((l != 0) && (ret == 0))
	{
		inverse = getBlock(t, 0, m * p, N + 2 * m, m, F);
		block = getBlock(0, k, 1, N, m, B);
		memcpy(inverse, block, l * sizeof(double));
		for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
		{
			block = getBlock(0, t + q * p, 1, N, m, B);
			temp = getBlock(t + q * p, k, N, N, m, A);
			matrixMultSubFrag(temp, inverse, block, m, l, 1);
		}
	}
	
	//debugDelete
	printf("B :\n");
	display(N, 1, m, B, N);


	return ret;
}


static inline int upStep_threads(double* A, double* B, double* F, int N, int m, int p, int t, int s, pthread_barrier_t* barrier)
{
	double *block, *inverse, *temp;

	pthread_barrier_wait(barrier);
	inverse = getBlock(t, 0, m * p, N + 2 * m, m, F);
	block = getBlock(0, s, 1, N, m, B);
	memcpy(inverse, block, m * sizeof(double));
	for(int q = 0; ((s - t - 1) >= 0) && (q <= (s - t - 1) / p); q++)
	{
		block = getBlock(0, t + q * p, 1, N, m, B);
		temp = getBlock(t + q * p, s, N, N, m, A);
		matrixMultSubFrag(temp, inverse, block, m, m, 1);
	}
	
	return 0;
}













static inline unsigned int solve_threads(int N, int m, int t, int p, double mtrNorm, double* A, double* B, double* F, pthread_barrier_t* barrier)
{
	int k = N / m;
	int err = 0;


	//sync
	pthread_barrier_wait(barrier);

	//Down
	for(int s = 0; s < k; s++)
	{
		//sync
		err = findCardinal_threads(p, t, A, B, N, m, mtrNorm, s, F, barrier);

		//debug
		//#ifdef __DEBUG__
		//	pthread_barrier_wait(barrier);
		//	if(t == 0)
		//	{
		//		printf("OOOOOOOOOOOOOOOOOOOOO:\n");
		//		display(N, N, m, A, N);
		//		display(N, 1, m, B, N);
		//	}
		//	pthread_barrier_wait(barrier);
		//#endif

		if(err < 0)
		{
			return (10 + s);
		}

		//sync
		multRow_threads(p, t, A, B, N, m, s, F, barrier);

		//debug
		//#ifdef __DEBUG__
		//	pthread_barrier_wait(barrier);
		//	if(t == 0)
		//	{
		//		printf("AAAAAAAAAAAAAAAAAAAAA:\n");
		//		display(N, N, m, A, N);
		//		display(N, 1, m, B, N);
		//	}
		//	pthread_barrier_wait(barrier);
		//#endif

		downStep_threads(A, B, F, N, m, p, t, s);

		//debug
		//#ifdef __DEBUG__
		//	pthread_barrier_wait(barrier);
		//	if(t == 0)
		//	{
		//		printf("BBBBBBBBBBBBBBBBBBBBB:\n");
		//		display(N, N, m, A, N);
		//		display(N, 1, m, B, N);
		//	}
		//	pthread_barrier_wait(barrier);
		//#endif

	}
	//sync
	
	
	#ifdef __DEBUG__
	if (t == p - 1)
	{
		printf("B almost :\n");
		display(N, 1, m, B, N);
	}
	#endif
	
	err = lastDownFirstUp_threads(A, B, N, m, p, t, mtrNorm, F, barrier);

	if(err < 0)
	{
		return (10 + k);
	}

	//debug
	//#ifdef __DEBUG__
	//	pthread_barrier_wait(barrier);
	//	if(t == 0)
	//	{
	//		display(N, N, m, A, N);
	//		display(N, 1, m, B, N);
	//		printf("////////////////////\n");
	//	}
	//	pthread_barrier_wait(barrier);
	//#endif

	//Up
	for(int s = k - 1; s > 0; s--)
	{
		//sync
		upStep_threads(A, B, F, N, m, p, t, s, barrier);

		//debug
		//#ifdef __DEBUG__
		//	pthread_barrier_wait(barrier);
		//	if(t == 0)
		//	{
		//		display(N, N, m, A, N);
		//		display(N, 1, m, B, N);
		//	}
		//	pthread_barrier_wait(barrier);
		//#endif
	}

	return 0;
}



void* thread_do(thread_args* args)
{
	unsigned int err = 0;

	int N = args->N;
	int m = args->m;
	int t =args->t;
	int p =args->p;
	double* A = args->A;
	double* B = args->B;
	double* X = args->X;
	double* F = args->F;
	double b;
	pthread_barrier_t* barrier = args->barrier;

	static double mtrNorm;


	//init
	err = init_thread(args->p, args->t, args->N, args->m, args->A, args->B, args->F, args->s, args->filename, args->barrier);
	if(err == 0)
	{
		if(args->t == 0)
		{
			display(args->N, args->N, args->m, args->A, args->r);
			display(args->N, 1, args->m, args->B, args->r);
			//display(args->N, 1, args->m, args->B, args->N);
			mtrNorm = norm(A, N, m, F);
			// printf("norm = %e\n", mtrNorm);
		}
		args->ret = 0;
		pthread_barrier_wait(args->barrier);
		args->thread_time = thread_time();
		args->full_time = full_time();
		args->ret = solve_threads(N, m, t, p, mtrNorm, A, B, F, barrier);
		args->thread_time = thread_time() - args->thread_time;
		args->full_time = full_time() - args->full_time;
		if(args->ret == 0)
		{
			reduce_max(&(args->full_time), barrier);
			if(t == 0)
				memcpy(X, B, (blockSizeD(1, m) * (N / m) + blockSizeD(1, (N - m * (N / m)))) * sizeof(double));
			pthread_barrier_wait(barrier);
			err = init_thread(p, t, N, m, A, B, F, args->s, args->filename, barrier);
			if(t == 0)
			{
				display(N, 1, m, X, args->r);
			}
			if((N <= 50) || (p > 1))
			{
				calculateResidual(A, B, X, F, N, m, p, t, &(args->res), &b);
				reduce_res(&(args->res), &b, barrier);
			}
			else
				args->res = -1.;
		}
	}
	else
	{
		args->ret = err;
	}
	//sync
	pthread_barrier_wait(args->barrier);
	return nullptr;
}

void* thread_func(void* vargs)
{
	int nprocs;
	cpu_set_t cpu;

	//cpuaffinity
	CPU_ZERO(&cpu);
	nprocs = get_nprocs();
	CPU_SET(nprocs - (((thread_args*)vargs)->t % nprocs), &cpu);
	sched_setaffinity(pthread_self(), sizeof(cpu), &cpu);

	thread_do((thread_args*)vargs);

	return nullptr;
}












//non3mult
// static int matrixMultFrag(double* A, double* B, double* C, int p, int q, int r)
// {
// 	int prem = p % 3, rrem = r % 3;

// 	for(int i = 0; i < p - prem; i += 3)
// 	{
// 		for(int j = 0; j < r - rrem; j += 3)
// 		{
// 			C[i * r + j] = 0.;
// 			C[i * r + j + 1] = 0.;
// 			C[i * r + j + 2] = 0.;
// 			C[(i + 1) * r + j] = 0.;
// 			C[(i + 1) * r + j + 1] = 0.;
// 			C[(i + 1) * r + j + 2] = 0.;
// 			C[(i + 2) * r + j] = 0.;
// 			C[(i + 2) * r + j + 1] = 0.;
// 			C[(i + 2) * r + j + 2] = 0.;
// 			for(int k = 0; k < q; k++)
// 			{
// 				C[i * r + j] += A[i * q + k] * B[k * r + j];
// 				C[i * r + j + 1] += A[i * q + k] * B[k * r + j + 1];
// 				C[i * r + j + 2]+= A[i * q + k] * B[k * r + j + 2];
// 				C[(i + 1) * r + j] += A[(i + 1) * q + k] * B[k * r + j];
// 				C[(i + 1) * r + j + 1] += A[(i + 1) * q + k] * B[k * r + j + 1];
// 				C[(i + 1) * r + j + 2] += A[(i + 1) * q + k] * B[k * r + j + 2];
// 				C[(i + 2) * r + j] += A[(i + 2) * q + k] * B[k * r + j];
// 				C[(i + 2) * r + j + 1] += A[(i + 2) * q + k] * B[k * r + j + 1];
// 				C[(i + 2) * r + j + 2] += A[(i + 2) * q + k] * B[k * r + j + 2];
// 			}
// 		}
// 		for(int j = r - rrem; j < r; j++)
// 		{
// 			C[i * r + j] = 0.;
// 			C[(i + 1) * r + j] = 0.;
// 			C[(i + 2) * r + j] = 0.;
// 			for(int k = 0; k < q; k++)
// 			{
// 				C[i * r + j] += A[i * q + k] * B[k * r + j];
// 				C[(i + 1) * r + j] += A[(i + 1) * q + k] * B[k * r + j];
// 				C[(i + 2) * r + j] += A[(i + 2) * q + k] * B[k * r + j];
// 			}
// 		}
// 	}
// 	for(int i = p - prem; i < p; i++)
// 	{
// 		for(int j = 0; j < r - rrem; j += 3)
// 		{
// 			C[i * r + j] = 0.;
// 			C[i * r + j + 1] = 0.;
// 			C[i * r + j + 2] = 0.;
// 			for(int k = 0; k < q; k++)
// 			{
// 				C[i * r + j] += A[i * q + k] * B[k * r + j];
// 				C[i * r + j + 1] += A[i * q + k] * B[k * r + j + 1];
// 				C[i * r + j + 2]+= A[i * q + k] * B[k * r + j + 2];
// 			}
// 		}
// 		for(int j = r - rrem; j < r; j++)
// 		{
// 			C[i * r + j] = 0.;
// 			for(int k = 0; k < q; k++)
// 			{
// 				C[i * r + j] += A[i * q + k] * B[k * r + j];
// 			}
// 		}
// 	}
// 	return 0;
// }


// static inline int matrixAdd(double* A, double* B, int p, int q, int sign)
// {
// 	if(sign < 0)
// 	{
// 		for(int i = 0; i < p; i++)
// 			for(int j = 0; j < q; j++)
// 				A[i * q + j] -= B[i * q + j];
// 	}
// 	else if(sign > 0)
// 	{
// 		for(int i = 0; i < p; i++)
// 			for(int j = 0; j < q; j++)
// 				A[i * q + j] += B[i * q + j];
// 	}
// 	else
// 	{
// 		for(int i = 0; i < p * q; i++)
// 			A[i] = -1.0 * A[i];
// 	}
// 	return 0;
// }
