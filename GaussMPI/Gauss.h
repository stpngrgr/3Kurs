#ifndef GAUSS_H_
#define GAUSS_H_

// #ifdef __DEBUG__
// #define __PROC_PRINT__ 0
// #define printproc (__PROC_PRINT__ < p ? __PROC_PRINT__ : 0)
// #else
#define __PROC_PRINT__ 0
#define printproc 0
// #endif


#include <pthread.h>
#include <math.h>
#include "mpi.h"






int init(double *A, double *B, double *F, int N, int m, int p, int t, int s, char *filename);


int solve_MPI(double *A, double *F, double *B, int N, int m, int p, int t, double cof);










static inline int blockSizeD(int rows, int cols)
{
	return rows * cols;
}

static inline int blockSizeB(int rows, int cols)
{
	return rows * cols * sizeof(double);
}


//static inline double* getBlock(int i, int j, int l, int n, int m, double* A)
//{
//	int kw = n / m, lw = n - kw * m;
//	int kh = l / m, lh = l - kh * m;
//	int fullsize = blockSizeD(m, m), fragsizew = blockSizeD(m, lw), fragsizeh = blockSizeD(m, lh);
//
//
//	return (A + i * (fullsize * kw + fragsizew) + j * (i < kh ? fullsize : fragsizeh));
//}



static inline double* getBlock(int i, int j, int l, int n, int m, double* A)
{
	int kw = n / m, lw = n - kw * m;
	int kh = l / m, lh = l - kh * m;
	int fullsize = blockSizeD(m, m), fragsizew = blockSizeD(m, lw), fragsizeh = blockSizeD(lh, m);


	return (A + j * (fullsize * kh + fragsizeh) + i * (j < kw ? fullsize : fragsizew));
}




static inline int awidth(int N, int m, int p)
{
	int k = N / m, l = N - k * m;

	return (k / p) + ((k % p != 0) || (l != 0) ? 1 : 0);
}


static inline int asizempi(int N, int m, int p)
{
	int k = N / m, l = N - k * m;
	int size = (blockSizeD(m, m) * k + blockSizeD(l, m)) * (k / p);

	if ((k % p != 0) || (l != 0))
	{
		size += (blockSizeD(m, m) * k + blockSizeD(l, m));
	}
	return size;
}


static inline int procawidth(int N, int m, int p, int t)
{
	int k = N / m, l = N - k * m;
	int ret = (k / p);

	if ((t < (k % p)) || ((t == (p - 1)) && (l != 0)))
		ret++;

	return ret;
}


static inline int procawidthfullonly(int N, int m, int p, int t)
{
	int k = N / m;
	// int l = N - k * m;
	int ret = (k / p);

	if (t < (k % p))
		ret++;

	return ret;
}


static inline int bsize(int N, int m)
{
	int k = N / m, l = N - k * m;

	return blockSizeD(m, 1) * k + blockSizeD(l, 1);
}


static inline int fsize(int N, int m)
{
	int k = N / m, l = N - k * m;

	return blockSizeD(m, m) * (k + 3) + blockSizeD(m, l) + 6;
}


static inline int colsize(int N, int m)
{
	int k = N / m, l = N - k * m;

	return blockSizeD(m, m) * k + blockSizeD(m, l);
}

static inline int fragcolsize(int N, int m)
{
	int k = N / m, l = N - k * m;

	return blockSizeD(m, m) * k + blockSizeD(m, l);
}




static inline int getColumnOwner(int N, int m, int p, int j)
{
	int k = N / m;

	if (j != k)
		return j % p;
	else
		return p - 1;
}


static inline int getColumnRank(int N, int m, int p, int j)
{
	int k = N / m;

	if (j != k)
		return j / p;
	else
		return procawidthfullonly(N, m, p, p - 1);
}




// static inline int L2GLinear(int L, int N, int m, int p, int t, int* Gblocki, int* Gblockj, int* Ginneri, int* Ginnerj)
// {
// 	int k = N / m, l = N - m * k;
// 	int colsize = k * blockSizeD(m, m) + blockSizeD(m, l);
// 	int inner;
// 	*Gblockj = L / colsize;
// 	if((t != (p - 1)) || (l == 0) || (lk != (k / p)))
// 	{
// 		*Gblocki = (L % colsize) / blockSizeD(m, m);
// 		inner = (L % colsize) - (*Gblocki) * blockSizeD(m, m);
// 		*Ginneri = inner / m;
// 		*Ginnerj = inner % m;
// 	}
// 	else
// 	{
// 		*Gblocki = (L % colsize) / blockSizeD(m, l);
// 		inner = (L % colsize) - (*Gblocki) * blockSizeD(m, l);
// 		*Ginneri = inner / l;
// 		*Ginnerj = inner % l;
// 	}
// 	return 0;
// }


static inline int L2GFull(int Lblocki, int Lblockj, int Linneri, int Linnerj, int N, int m, int p, int t, int* Gblocki, int* Gblockj, int* Ginneri, int* Ginnerj)
{
	int k = N / m, l = N - m * k;
	// int colsize = k * blockSizeD(m, m) + blockSizeD(m, l);
	
	*Gblocki = Lblocki;
	*Ginneri = Linneri;
	*Ginnerj = Linnerj;
	if ((t != (p - 1)) || (l == 0) || (Lblockj != procawidth(N, m, p, t)))
	{
		*Gblockj = Lblockj * p + t;
	}
	else
	{
		*Gblockj = k;
	}
	return 0;
}

static inline int L2GBlock(int Lblocki, int Lblockj, int N, int m, int p, int t, int* Gblocki, int* Gblockj)
{
	int k = N / m, l = N - m * k;
	// int colsize = k * blockSizeD(m, m) + blockSizeD(m, l);
	*Gblocki = Lblocki;
	if ((t != (p - 1)) || (l == 0) || (Lblockj != getColumnRank(N, m, p, k)))
	{
		*Gblockj = Lblockj * p + t;
	}
	else
	{
		*Gblockj = k;
	}
	return 0;
}

static inline int L2GBJ(int L, int N, int m, int p, int t)
{
	int k = N / m, l = N - m * k;
	int ret;
	if ((t != (p - 1)) || (l == 0) || (L != procawidth(N, m, p, t)))
	{
		ret = L * p + t;
	}
	else
	{
		ret = k;
	}
	return ret;
}

static inline int L2GBI(int L, int N, int m, int p, int t)
{
	(void)(N);
	(void)(m);
	(void)(p);
	(void)(t);
	return L;
}

//static inline int G2L(int G, int N, int m, int p)
//{
//	int L;
//	int k = N / m, l = N - m * k;
//	int lk = lj / (m * N);
//	if((t != (p - 1)) || (l == 0) || (lk != (k / p)))
//		L = lk * (m * N * (p + t)) + (lj % (m * N));
//	else
//		L = N * m * k + (lj % (m * N));
//	return L;
//}



















static inline int displayA_MPI(double* A, double* F, int N, int m, int p, int t, int r)
{
	int k = N / m, l = N - k * m;
	int printsize = (r < N ? r : N);
	int blockwidth, printwidth;
	int printblocksize = (printsize / m) + (printsize % m ? 1 : 0);
	int owner;
	double *block;
	MPI_Comm world = MPI_COMM_WORLD;
	MPI_Status status;


	if (t == printproc)
		printf("\n\n");
	for (int i = 0; i < printsize; i++)
	{
		for (int blockj = 0; blockj < printblocksize; blockj++)
		{
			owner = getColumnOwner(N, m, p, blockj);
			if (t == printproc)
			{
				if (t == owner)
					block = getBlock(i / m, getColumnRank(N, m, p, blockj), N, awidth(N, m, p) * m, m, A);
				else
				{
					MPI_Recv(F, m * m, MPI_DOUBLE, owner, 0, world, &status);
					block = F;
				}
				blockwidth = (blockj != k ? m : l);
				printwidth = (blockj * m + blockwidth < printsize ? blockwidth : printsize - blockj * m);
				for (int j = 0; j < printwidth; j++)
				{
					printf(" %10.3e", block[(i % m) * blockwidth + j]);
				}
			}
			else if (t == owner)
			{
				block = getBlock(i / m, getColumnRank(N, m, p, blockj), N, awidth(N, m, p) * m, m, A);
				MPI_Send(block, m * m, MPI_DOUBLE, printproc, 0, world);
			}
		}
		if (t == printproc)
			printf("\n");
	}
	if (t == printproc)
		printf("\n\n");
	return 0;
}




static inline int display(double* A, int l, int n, int m, int r)
{
	int kw = n / m;
	int lw = n - kw * m;
	// int kh = l / m;
	// int lh = l - kh * m;
	int h = (l > r ? r : l), w = (n > r ? r : n);
	double p, *block;

	printf("\n\n");
	for(int i = 0; i < h; i++)
	{
 		for(int j = 0; j < w; j++)
 		{
 			block = getBlock((i / m), (j / m), l, n, m, A);
 			p = block[(i % m) * ((j - kw * m) < 0 ? m : lw) + (j % m)];
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




#endif /* GAUSS_H_ */
