#include "Gauss.h"

//#include <sys/sysinfo.h>
//#include <sched.h>
//#include <sys/resource.h>
//#include <sys/time.h>
#include <stdio.h>
#include <pthread.h>
#include <string.h>

/*matrixNorm_threads(p, t, F, N, barrier);
findCardinal_threads(p, t, A, B, N, m, s, ret, F);
multRow_threads(p, t, N, m, s, F);
downStep_threads(A, B, F, N, m, p, t, s);
lastDownFirstUp_threads(A, B, N, m, p, t, F, ret);
upStep_threads(A, B, F, N, m, p, t, s);
init_threads(A, N, m, p, t, sf, filename, barrier);*/


int blockSizeD(int rows, int cols)
{
	return rows * cols;
}

int blockSizeB(int rows, int cols)
{
	return blockSizeD(rows, cols) * sizeof(double);
}

double* getBlock(int i, int j, int l, int n, int m, double* A)
{
	int kw = n / m, lw = n - kw * m;
	int kh = l / m, lh = l - kh * m;
	int fullsize = blockSizeD(m, m), fragsizew = blockSizeD(m, lw), fragsizeh = blockSizeD(m, lh);


	return (A + i * (fullsize * kw + fragsizew) + j * (i < kh ? fullsize : fragsizeh));
}

void dcpy(double* dest, double* from, int size)
{
	for(int i = 0; i < size; i++)
	{
		dest[i] = from[i];
	}
}

void mcpy(void* dest, void* from, int size)
{
	char *d = (char*)dest, *f = (char*)from;


	for(int i = 0; i < size; i++)
	{
		d[i] = f[i];
	}
}

int args_init(thread_args* args, int p, int N, int m, double* A, double* B, double* X, double *F, int s, int r, char* filename, pthread_barrier_t* barrier)
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

void reduce_err(int* ret, pthread_barrier_t* barrier)
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

void reduce_erru(unsigned int* ret, pthread_barrier_t* barrier)
{
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static unsigned int global = 0;
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

void reduce_max(double* s, pthread_barrier_t* barrier)
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

void reduce_min_di(double* d, int* i, pthread_barrier_t* barrier)
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




int findCardinal_threads(int p, int t, double* A, double* B, int N, int m, double mtrNorm, int s, double* F, pthread_barrier_t* barrier)
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
		mcpy(Fcopy, Ablock, m * m * sizeof(double));
		err = matrixInverse(Fcopy, Ftemp, m, mtrNorm);
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
				if(((k - s - 1) - ((k - s - 1) / 2)) > 0)
				{
					Ablock = getBlock(s, s + 1, N, N, m, A);
					Ftemp = getBlock(mei, s + 1, N, N, m, A);
					Fcopy = getBlock(t, 1, m * p, N + 2 * m, m, F);
					memcpy(Fcopy, Ftemp, (blockSizeD(m, m) * ((k - s - 1) - ((k - s - 1) / 2))) * sizeof(double));
					memcpy(Ftemp, Ablock, (blockSizeD(m, m) * ((k - s - 1) - ((k - s - 1) / 2))) * sizeof(double));
					memcpy(Ablock, Fcopy, (blockSizeD(m, m) * ((k - s - 1) - ((k - s - 1) / 2))) * sizeof(double));
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
			if((s != mei) && ((l != 0) || (((k - s - 1) / 2) > 0)))
			{
				Ablock = getBlock(s, s + 1 + ((k - s - 1) - ((k - s - 1) / 2)), N, N, m, A);
				Ftemp = getBlock(mei, s + 1 + ((k - s - 1) - ((k - s - 1) / 2)), N, N, m, A);
				Fcopy = getBlock(t, 1, m * p, N + 2 * m, m, F);
				memcpy(Fcopy, Ftemp, (blockSizeD(m, m) * ((k - s - 1) / 2) + blockSizeD(m, l)) * sizeof(double));
				memcpy(Ftemp, Ablock, (blockSizeD(m, m) * ((k - s - 1) / 2) + blockSizeD(m, l)) * sizeof(double));
				memcpy(Ablock, Fcopy, (blockSizeD(m, m) * ((k - s - 1) / 2) + blockSizeD(m, l)) * sizeof(double));
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
		err = 0;
	}
	else
	{
		err = -1;
	}
	return err;
}

int multRow_threads(int p, int t, double* A, double* B, int N, int m, int s, double* F, pthread_barrier_t* barrier)
{
	int k = N / m, l = N - k * m;
	double *block, *inverse, *temp;


	inverse = getBlock(t, k, m * p, N + 2 * m, m, F);

	//debug
	// if(t == 0)
	// 	display(m, m, m, inverse, m);

	temp = getBlock(t, k + 1, m * p, N + 2 * m, m, F);
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
	if(t == (s % p))
	{
		block = getBlock(0, s, 1, N, m, B);
		matrixMultRewriteRight(inverse, block, temp, m, 1);
	}

	//sync
	pthread_barrier_wait(barrier);

	for(int ti = 0; ti < p; ti++)
	{
		if(ti != t)
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
	}
	block = getBlock(t, k, m * p, N + 2 * m, m, F);
	temp = getBlock(0, s, 1, N, m, B);
	memcpy(block, temp, m * sizeof(double));
	if(t == (s % p))
	{
		if(s < (k - 1))
		{
			block = getBlock(s, s + 1, N, N, m, A);
			temp = getBlock(t, s + 1, m * p, N + 2 * m, m, F);
			memcpy(block, temp, blockSizeD(m, m) * (k - s - 1) * sizeof(double));
		}
		block = getBlock(s, k, N, N, m, A);
		temp = getBlock(t, k + 2, m * p, N + 2 * m, m, F);
		memcpy(block, temp, blockSizeD(m, l) * sizeof(double));
	}

	return 0;
}

int downStep_threads(double* A, double* B, double* F, int N, int m, int p, int t, int s)
{
	int k = N / m, l = N - k * m;
	double *AblockMult, *AblockSub, *Fblock;


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
			matrixMultSubFrag(AblockMult, Fblock, AblockSub, l, m, l);
		}
		Fblock = getBlock(t, k, m * p, N + 2 * m, m, F);
		AblockSub = getBlock(0, k, 1, N, m, B);
		matrixMultSubFrag(AblockMult, Fblock, AblockSub, l, m, 1);
	}

	return 0;
}

int lastDownFirstUp_threads(double* A, double* B, int N, int m, int p, int t, double norm, double* F, pthread_barrier_t* barrier)
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

			ret = matrixInverse(block, inverse, l, norm);
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
		dcpy(inverse, block, l);
		for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
		{
			block = getBlock(0, t + q * p, 1, N, m, B);
			temp = getBlock(t + q * p, k, N, N, m, A);
			matrixMultSubFrag(temp, inverse, block, m, l, 1);
		}
	}

	return ret;
}


int upStep_threads(double* A, double* B, double* F, int N, int m, int p, int t, int s, pthread_barrier_t* barrier)
{
	double *block, *inverse, *temp;

	pthread_barrier_wait(barrier);
	inverse = getBlock(t, 0, m * p, N + 2 * m, m, F);
	block = getBlock(0, s, 1, N, m, B);
	dcpy(inverse, block, m);
	for(int q = 0; ((s - t - 1) >= 0) && (q <= (s - t - 1) / p); q++)
	{
		block = getBlock(0, t + q * p, 1, N, m, B);
		temp = getBlock(t + q * p, s, N, N, m, A);
		matrixMultSubFrag(temp, inverse, block, m, m, 1);
	}
	
	return 0;
}


/*template1

int k = N / m, l = N - k * m;
double *block;

for(int q = (s / p) + ((s % p) > t ? 1 : 0); ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
{

}
if((t == (p - 1)) && (l != 0))
{

}

*/


//int matrixNorm_threads(int p, int t, int N, int m, double* A, double* F, pthread_barrier_t* barrier)
//{
//	int k = N / m, l = N - k * m;
//	int ret = 0, barret;
//	double* Ablock, Fblock;
//
//	static double res;
//
//	//local
//	for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p); q++)
//	{
//		for(int b = 0; b < k; b++)
//		{
//			Fblock = getBlock(t, b, N, m, F);
//			Ablock = getBlock(t + q * p, b, N, m, A);
//			for(int j = 0; j < m; j++)
//			{
//				Fblock[j] = 0.;
//				for(int i = 0; i < m; i++)
//				{
//					Fblock[j] += Ablock[i * m + j];
//				}
//			}
//		}
//		if(l != 0)
//		{
//			Fblock = getBlock(t, k, N, m, F);
//			Ablock = getBlock(t + q * p, k, N, m, A);
//			for(int j = 0; j < m; j++)
//			{
//				Fblock[j] = 0.;
//				for(int i = 0; i < l; i++)
//				{
//					Fblock[j] += Ablock[i * m + j];
//				}
//			}
//
//		}
//	}
//	if((t == (p - 1)) && (l != 0))
//	{
//		for(int b = 0; b < k; b++)
//		{
//			Fblock = getBlock(p - 1, b, N, m, F);
//			Ablock = getBlock(k, b, N, m, A);
//			for(int j = 0; j < m; j++)
//			{
//				Fblock[j] = 0.;
//				for(int i = 0; i < m; i++)
//				{
//					Fblock[j] += Ablock[i * m + j];
//				}
//			}
//		}
//		Fblock = getBlock(p - 1, k, N, m, F);
//		Ablock = getBlock(k, k, N, m, A);
//		for(int j = 0; j < l; j++)
//		{
//			Fblock[j] = 0.;
//			for(int i = 0; i < l; i++)
//			{
//				Fblock[j] += Ablock[i * m + j];
//			}
//		}
//	}
//
//	//sync
//	barret = pthread_barrier_wait(barrier);
//	if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
//	{
//		for(int b = 0; b < k; b++)
//		{
//			for(int i = 0; i < m; i++)
//			{
//				for(int tt = 0; tt < p; tt++)
//				{
//
//				}
//			}
//		}
//		if(l != 0)
//		{
//
//		}
//	}
//}
