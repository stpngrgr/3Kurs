#include <stdio.h>
#include <math.h>
#include <string.h>
#include "Gauss.h"

//#include <sys/sysinfo.h>
//#include <sched.h>
//#include <sys/resource.h>
#include <pthread.h>




int init_formula(double *A, double *B, int N, int m, int p, int t, int s);
int init_file(double* A, double* B, double* buffer, int N, int m, int p, int t, char* filename);
int readA(double* A, int N, int m, FILE* f, int step);
int bufcpy(double* A, double* buffer, int N, int m, int p, int t, int i);
int updB(double* B, double* buffer, int N, int m, int i);






int init(double *A, double *B, double *F, int N, int m, int p, int t, int s, char *filename)
{
	memset(A, 0, asizempi(N, m, p) * sizeof(double));
	memset(B, 0, bsize(N, m) * sizeof(double));
	memset(F, 0, (fsize(N, m) + bsize(N, m)) * sizeof(double));
	if (s == 0)
		return init_file(A, B, F, N, m, p, t, filename);
	else
		return init_formula(A, B, N, m, p, t, s);
}



int init_formula(double *A, double *B, int N, int m, int p, int t, int s)
{
	int k = N / m, l = N - k * m;
	int bcols, brows = k;
	double *block, *bblock = nullptr;
	int blockrownum, blockcolnum;
	int Gblocki, Gblockj;
	int err = 0;

	bcols = procawidth(N, m, p, t);
	if (l != 0)
	{
		brows++;
	}
	
	switch(s)
	{
	case 1 : //n - max{i, j} + 1
		for (int blj = 0; blj < bcols; blj++)
		{
			for (int bli = 0; bli < brows; bli++)
			{
				L2GBlock(bli, blj, N, m, p, t, &Gblocki, &Gblockj);
				block = getBlock(bli, blj, N, awidth(N, m, p) * m, m, A);
				blockcolnum = (Gblockj != k ? m : l);
				blockrownum = (Gblocki != k ? m : l);
				for (int i = 0; i < blockrownum; i++)
				{
					for (int j = 0; j < blockcolnum; j++)
					{
						block[i * blockcolnum + j] = N - (Gblocki * m + i > Gblockj * m + j ? Gblocki * m + i : Gblockj * m + j);
					}
				}
			}
		}
		for (int bli = 0; bli < brows; bli++)
		{
			bblock = getBlock(bli, 0, N, 1, m, B);
			blockrownum = (bli != k ? m : l);
			for (int i = 0; i < blockrownum; i++)
				for (int j = 0; j < N; j += 2)
				{
					bblock[i] += N - (bli * m + i > j ? bli * m + i : j);
				}
		}
		break;
	case 2 : //max{i, j}
		for (int blj = 0; blj < bcols; blj++)
		{
			for (int bli = 0; bli < brows; bli++)
			{
				L2GBlock(bli, blj, N, m, p, t, &Gblocki, &Gblockj);
				block = getBlock(bli, blj, N, awidth(N, m, p) * m, m, A);
				blockcolnum = (Gblockj != k ? m : l);
				blockrownum = (Gblocki != k ? m : l);
				for (int i = 0; i < blockrownum; i++)
				{
					for (int j = 0; j < blockcolnum; j++)
					{
						block[i * blockcolnum + j] = (Gblocki * m + i > Gblockj * m + j ? Gblocki * m + i : Gblockj * m + j) + 1;
					}
				}
			}
		}
		for (int bli = 0; bli < brows; bli++)
		{
			bblock = getBlock(bli, 0, N, 1, m, B);
			blockrownum = (bli != k ? m : l);
			for (int i = 0; i < blockrownum; i++)
				for (int j = 0; j < N; j += 2)
				{
					bblock[i] += (bli * m + i > j ? bli * m + i : j) + 1;
				}
		}
		break;
	case 3 : //|i - j|
		for (int blj = 0; blj < bcols; blj++)
		{
			for (int bli = 0; bli < brows; bli++)
			{
				L2GBlock(bli, blj, N, m, p, t, &Gblocki, &Gblockj);
				block = getBlock(bli, blj, N, awidth(N, m, p) * m, m, A);
				blockcolnum = (Gblockj != k ? m : l);
				blockrownum = (Gblocki != k ? m : l);
				for (int i = 0; i < blockrownum; i++)
				{
					for (int j = 0; j < blockcolnum; j++)
					{
						block[i * blockcolnum + j] = abs((Gblocki * m + i) - (Gblockj * m + j));
					}
				}
			}
		}
		for (int bli = 0; bli < brows; bli++)
		{
			bblock = getBlock(bli, 0, N, 1, m, B);
			blockrownum = (bli != k ? m : l);
			for (int i = 0; i < blockrownum; i++)
				for (int j = 0; j < N; j += 2)
				{
					bblock[i] += abs((bli * m + i) - j);
				}
		}
		break;
	case 4 : //1/(i + j - 1)
		for (int blj = 0; blj < bcols; blj++)
		{
			for (int bli = 0; bli < brows; bli++)
			{
				L2GBlock(bli, blj, N, m, p, t, &Gblocki, &Gblockj);
				block = getBlock(bli, blj, N, awidth(N, m, p) * m, m, A);
				blockcolnum = (Gblockj != k ? m : l);
				blockrownum = (Gblocki != k ? m : l);
				for (int i = 0; i < blockrownum; i++)
				{
					for (int j = 0; j < blockcolnum; j++)
					{
						block[i * blockcolnum + j] = 1. / ((Gblocki * m + i) + (Gblockj * m + j) + 1);
					}
				}
			}
		}
		for (int bli = 0; bli < brows; bli++)
		{
			bblock = getBlock(bli, 0, N, 1, m, B);
			blockrownum = (bli != k ? m : l);
			for (int i = 0; i < blockrownum; i++)
				for (int j = 0; j < N; j += 2)
				{
					bblock[i] += 1. / ((bli * m + i) + j + 1);
				}
		}
		break;
	default :
		err = -1;
		break;
	}

	return err;
}




int init_file(double* A, double* B, double* buffer, int N, int m, int p, int t, char* filename)
{
	int k = N / m, l = N - k * m;
	MPI_Comm world = MPI_COMM_WORLD;
	int err = 0;
	FILE* f = nullptr;

	if (t == 0)
	{
		f = fopen(filename, "r");
		if (f == nullptr)
			err = -2;
	}
	MPI_Bcast(&err, 1, MPI_INT, 0, world);
	if (err)
		return err;
	for (int i = 0; i < k + (l != 0 ? 1 : 0); i++)
	{
		if (t == 0)
		{
			if (!err)
			{
				err = readA(buffer, N, m, f, i);
				//debug
				// {
				// 	printf("\n %d \n", err);
				// 	if (t == 0)
				// 		for (int o = 0; o < colsize(N, m); o++)
				// 		{
				// 			printf(" %10.3e", buffer[o]);
				// 		}
				// 	printf("\n");
				// }
			}
		}
		MPI_Bcast(buffer, colsize(N, m), MPI_DOUBLE, 0, world);
		bufcpy(A, buffer, N, m, p, t, i);
		if (t == (p - 1))
		{
			updB(B, buffer, N, m, i);
		}
	}
	MPI_Bcast(B, bsize(N, m), MPI_DOUBLE, p - 1, world);
	MPI_Bcast(&err, 1, MPI_INT, 0, world);
	return err;
}




int readA(double* A, int N, int m, FILE* f, int step)
{
	int k = N / m, l = N - k * m;
	int im = ((step != k)? m : l);
	int blockWidth = 0, ret = 0;
	double *block = A;

	for (int i = 0; i < im; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if ((j % m) == 0)
			{
				block = getBlock(j / m, 0, N, m, m, A);
				blockWidth = ((j / m) != k ? m : l);
			}
			ret = fscanf(f, "%lf", block + i * blockWidth + (j % m));
			if (ret != 1)
				return (!feof(f) ? -3 : -4);
		}
	}
	return 0;
}



int bufcpy(double* A, double* buffer, int N, int m, int p, int t, int i)
{
	int k = N / m, l = N - m * k;
	int bnum = (k / p) + (t < (k % p) ? 1 : 0);
	int AW = awidth(N, m, p);
	int blockh = (i < k ? m : l), blockw = m;
	double *Ablock = A, *Bufblock = buffer;

	for (int j = 0; j < bnum; j++)
	{
		Ablock = getBlock(i, j, N, AW * m, m, A);
		Bufblock = getBlock(j * p + t, 0, N, m, m, buffer);
		memcpy(Ablock, Bufblock, blockSizeB(blockh, blockw));
	}
	if ((t == (p - 1)) && (l != 0))
	{
		blockw = l;
		Ablock = getBlock(i, (k / p), N, AW * m, m, A);
		Bufblock = getBlock(k, 0, N, m, m, buffer);
		memcpy(Ablock, Bufblock, blockSizeB(blockh, blockw));
	}
	return 0;
}



int updB(double* B, double* buffer, int N, int m, int i)
{
	int k = N / m, l = N - k * m;
	int blockh = (i < k ? m : l);
	double *Bblock = getBlock(i, 0, N, 1, m, B), *Bufblock;

	for (int blockj = 0; blockj < k; blockj++)
	{
		Bufblock = getBlock(0, blockj, m, N, m, buffer);
		for (int inneri = 0; inneri < blockh; inneri++)
			for (int innerj = ((blockj * m) % 2 ? 1 : 0); innerj < m; innerj += 2)
				Bblock[inneri] += Bufblock[inneri * m + innerj];
	}
	if (l != 0)
	{
		Bufblock = getBlock(0, k, m, N, m, buffer);
		for (int inneri = 0; inneri < blockh; inneri++)
			for (int innerj = ((k * m) % 2 ? 1 : 0); innerj < l; innerj += 2)
				Bblock[inneri] += Bufblock[inneri * l + innerj];
	}
	return 0;
}

