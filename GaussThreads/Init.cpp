#include <stdio.h>
#include <math.h>
#include "Gauss.h"

//#include <sys/sysinfo.h>
//#include <sched.h>
//#include <sys/resource.h>
#include <pthread.h>



static inline void reduce_erru(unsigned int* ret, pthread_barrier_t* barrier)
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


unsigned int init_thread(int p, int t, int N, int m, double* A, double* B, double* F, int s, char* filename, pthread_barrier_t* barrier)
{
	int k = N / m, l = N - k * m;
	int barret;
	unsigned int ret = 0;
	int ii, jj;
	double *block, *block1;
	FILE* f;

	//sync
	pthread_barrier_wait(barrier);
	if(A == nullptr)
		ret = 1;
	//sync
	reduce_erru(&ret, barrier);
	if(ret == 0)
	{
		switch(s)
		{
		case 0:

			for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p)  ; q++)
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(t + q * p, j, N, N, m, A);
					for(int i = 0; i < blockSizeD(m, m); i++)
						block[i] = 0.;
				}
				if(l != 0)
				{
					block = getBlock(t + q * p, k, N, N, m, A);
					for(int i = 0; i < blockSizeD(m, l); i++)
						block[i] = 0.;
				}
			}
			if((t == (p - 1)) && (l != 0))
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(k, j, N, N, m, A);
					for(int i = 0; i < m * l; i++)
						block[i] = 0.;
				}
				block = getBlock(k, k, N, N, m, A);
				for(int i = 0; i < l * l; i++)
					block[i] = 0.;
			}
			//sync
			barret = pthread_barrier_wait(barrier);
			if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
			{
				if((f = fopen(filename, "r")) == nullptr)
				{
					ret = 2;
				}
				if(ret == 0)
				{
					for(int i = 0; i < N; i++)
					{
						for(int j = 0; j < N; j++)
						{
							block = getBlock(i / m, j / m, N, N, m, A);
							if(1 != fscanf( f, "%lf", block + (i % m) * ((j / m) < k ? m : l) + (j % m)))
							{
								ret = (!feof(f) ? 3 : 4);
								break;
							}
						}
						if(ret != 0)
							break;
					}
					fclose(f);
				}
			}
			break;

		case 1:

			for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p)  ; q++)
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(t + q * p, j, N, N, m, A);
					for(int i = 0; i < m * m; i++)
					{
						ii = (t + q * p) * m + (i / m);
						jj = j * m + (i % m);
						block[i] = N - (ii > jj ? ii : jj);
					}
					for(int i = m * m; i < blockSizeD(m, m); i++)
					{
						block[i] = 0.;
					}
				}
				if(l != 0)
				{
					block = getBlock(t + q * p, k, N, N, m, A);
					for(int i = 0; i < m * l; i++)
					{
						ii = (t + q * p) * m + (i / l);
						jj = k * m + (i % l);
						block[i] = N - (ii > jj ? ii : jj);
					}
					for(int i = m * l; i < blockSizeD(m, l); i++)
					{
						block[i] = 0.;
					}
				}
			}
			if((t == (p - 1)) && (l != 0))
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(k, j, N, N, m, A);
					for(int i = 0; i < m * l; i++)
					{
						ii = k * m + (i / m);
						jj = j * m + (i % m);
						block[i] = N - (ii > jj ? ii : jj);
					}
					for(int i = m * l; i < blockSizeD(m, l); i++)
					{
						block[i] = 0.;
					}
				}
				block = getBlock(k, k, N, N, m, A);
				for(int i = 0; i < l * l; i++)
				{
					ii = k * m + (i / l);
					jj = k * m + (i % l);
					block[i] = N - (ii > jj ? ii : jj);
				}
				for(int i = l * l; i < blockSizeD(l, l); i++)
				{
					block[i] = 0.;
				}
			}
			break;

		case 2:

			for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p)  ; q++)
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(t + q * p, j, N, N, m, A);
					for(int i = 0; i < m * m; i++)
					{
						ii = (t + q * p) * m + (i / m);
						jj = j * m + (i % m);
						block[i] = (ii > jj ? ii : jj) + 1;
					}
					for(int i = m * m; i < blockSizeD(m, m); i++)
					{
						block[i] = 0.;
					}
				}
				if(l != 0)
				{
					block = getBlock(t + q * p, k, N, N, m, A);
					for(int i = 0; i < m * l; i++)
					{
						ii = (t + q * p) * m + (i / l);
						jj = k * m + (i % l);
						block[i] = (ii > jj ? ii : jj) + 1;
					}
					for(int i = m * l; i < blockSizeD(m, l); i++)
					{
						block[i] = 0.;
					}
				}
			}
			if((t == (p - 1)) && (l != 0))
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(k, j, N, N, m, A);
					for(int i = 0; i < m * l; i++)
					{
						ii = k * m + (i / m);
						jj = j * m + (i % m);
						block[i] = (ii > jj ? ii : jj) + 1;
					}
					for(int i = m * l; i < blockSizeD(m, l); i++)
					{
						block[i] = 0.;
					}
				}
				block = getBlock(k, k, N, N, m, A);
				for(int i = 0; i < l * l; i++)
				{
					ii = k * m + (i / l);
					jj = k * m + (i % l);
					block[i] = (ii > jj ? ii : jj) + 1;
				}
				for(int i = l * l; i < blockSizeD(l, l); i++)
				{
					block[i] = 0.;
				}
			}
			break;

		case 3:


			for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p)  ; q++)
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(t + q * p, j, N, N, m, A);
					for(int i = 0; i < m * m; i++)
					{
						ii = (t + q * p) * m + (i / m);
						jj = j * m + (i % m);
						block[i] = abs(ii - jj);
					}
					for(int i = m * m; i < blockSizeD(m, m); i++)
					{
						block[i] = 0.;
					}
				}
				if(l != 0)
				{
					block = getBlock(t + q * p, k, N, N, m, A);
					for(int i = 0; i < m * l; i++)
					{
						ii = (t + q * p) * m + (i / l);
						jj = k * m + (i % l);
						block[i] = abs(ii - jj);
					}
					for(int i = m * l; i < blockSizeD(m, l); i++)
					{
						block[i] = 0.;
					}
				}
			}
			if((t == (p - 1)) && (l != 0))
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(k, j, N, N, m, A);
					for(int i = 0; i < m * l; i++)
					{
						ii = k * m + (i / m);
						jj = j * m + (i % m);
						block[i] = abs(ii - jj);
					}
					for(int i = m * l; i < blockSizeD(m, l); i++)
					{
						block[i] = 0.;
					}
				}
				block = getBlock(k, k, N, N, m, A);
				for(int i = 0; i < l * l; i++)
				{
					ii = k * m + (i / l);
					jj = k * m + (i % l);
					block[i] = abs(ii - jj);
				}
				for(int i = l * l; i < blockSizeD(l, l); i++)
				{
					block[i] = 0.;
				}
			}
			break;

		case 4:


			for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p)  ; q++)
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(t + q * p, j, N, N, m, A);
					for(int i = 0; i < m * m; i++)
					{
						ii = (t + q * p) * m + (i / m);
						jj = j * m + (i % m);
						block[i] = 1./(ii + jj + 1);
					}
					for(int i = m * m; i < blockSizeD(m, m); i++)
					{
						block[i] = 0.;
					}
				}
				if(l != 0)
				{
					block = getBlock(t + q * p, k, N, N, m, A);
					for(int i = 0; i < m * l; i++)
					{
						ii = (t + q * p) * m + (i / l);
						jj = k * m + (i % l);
						block[i] = 1./(ii + jj + 1);
					}
					for(int i = m * l; i < blockSizeD(m, l); i++)
					{
						block[i] = 0.;
					}
				}
			}
			if((t == (p - 1)) && (l != 0))
			{
				for(int j = 0; j < k; j++)
				{
					block = getBlock(k, j, N, N, m, A);
					for(int i = 0; i < m * l; i++)
					{
						ii = k * m + (i / m);
						jj = j * m + (i % m);
						block[i] = 1./(ii + jj + 1);
					}
					for(int i = m * l; i < blockSizeD(m, l); i++)
					{
						block[i] = 0.;
					}
				}
				block = getBlock(k, k, N, N, m, A);
				for(int i = 0; i < l * l; i++)
				{
					ii = k * m + (i / l);
					jj = k * m + (i % l);
					block[i] = 1./(ii + jj + 1);
				}
				for(int i = l * l; i < blockSizeD(l, l); i++)
				{
					block[i] = 0.;
				}
			}
			break;

		default:
			ret = 9;
		}

		//sync
		reduce_erru(&ret, barrier);

		if(ret == 0)
		{
			//B initialization
			for(int q = 0; ((k - t - 1) >= 0) && (q <= (k - t - 1) / p) ; q++)
			{
				block = getBlock(0, t + q * p, 1, N, m, B);
				for(int i = 0; i < blockSizeD(1, m); i++)
				{
					block[i] = 0.;
				}
				for(int jb = 0; jb < k; jb++)
				{
					block1 = getBlock(t + q * p, jb, N, N, m, A);
					for(int i = 0; i < m; i++)
					{
						for(int j = (((m * jb) % 2) == 0 ? 0 : 1); j < m; j += 2)
						{
							//printf("chto: %lf, %lf\n", block[i], block1[i * m + j]);
							block[i] += block1[i * m + j];
						}
					}
				}
				if(l != 0)
				{
					block1 = getBlock(t + q * p, k, N, N, m, A);
					for(int i = 0; i < m; i++)
					{
						for(int j = (((m * k) % 2) == 0 ? 0 : 1); j < l; j += 2)
						{
							block[i] += block1[i * l + j];
						}
					}
				}
			}
			if((t == (p - 1)) && (l != 0))
			{
				block = getBlock(0, k, 1, N, m, B);
				for(int i = 0; i < blockSizeD(1, l); i++)
				{
					block[i] = 0.;
				}
				for(int jb = 0; jb < k; jb++)
				{
					block1 = getBlock(k, jb, N, N, m, A);
					for(int i = 0; i < l; i++)
					{
						for(int j = (((m * jb) % 2) == 0 ? 0 : 1); j < m; j += 2)
						{
							block[i] += block1[i * m + j];
						}
					}
				}
				block1 = getBlock(k, k, N, N, m, A);
				for(int i = 0; i < l; i++)
				{
					for(int j = (((m * k) % 2) == 0 ? 0 : 1); j < l; j += 2)
					{
						block[i] += block1[i * l + j];
					}
				}
			}
			for(int jb = 0; jb < k + 2; jb++)
			{
				block = getBlock(t, jb, m * p, N + 2 * m, m, F);
				for(int i = 0; i < blockSizeD(m, m); i++)
				{
					block[i] = 0.;
				}
			}
			if(l != 0)
			{
				block = getBlock(t, k + 2, m * p, N + 2 * m, m, F);
				for(int i = 0; i < blockSizeD(m, l); i++)
				{
					block[i] = 0.;
				}
			}
		}
	}

	return ret;
}


