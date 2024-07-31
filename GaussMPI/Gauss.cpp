#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Gauss.h"


#include <stddef.h>
#include <fenv.h>
#include <unistd.h>


int main_do(int argc, char** argv, double **A, double **B, double **F, double **X, double **proctime);







static inline int matrixMultFrag(double* A, double* B, double* C, int p, int q, int r)
{
	int prem = p % 3, rrem = r % 3;

	for(int i = 0; i < p - prem; i += 3)
	{
		for(int j = 0; j < r - rrem; j += 3)
		{
			C[i * r + j] = 0.;
			C[i * r + j + 1] = 0.;
			C[i * r + j + 2] = 0.;
			C[(i + 1) * r + j] = 0.;
			C[(i + 1) * r + j + 1] = 0.;
			C[(i + 1) * r + j + 2] = 0.;
			C[(i + 2) * r + j] = 0.;
			C[(i + 2) * r + j + 1] = 0.;
			C[(i + 2) * r + j + 2] = 0.;
			for(int k = 0; k < q; k++)
			{
				C[i * r + j] += A[i * q + k] * B[k * r + j];
				C[i * r + j + 1] += A[i * q + k] * B[k * r + j + 1];
				C[i * r + j + 2]+= A[i * q + k] * B[k * r + j + 2];
				C[(i + 1) * r + j] += A[(i + 1) * q + k] * B[k * r + j];
				C[(i + 1) * r + j + 1] += A[(i + 1) * q + k] * B[k * r + j + 1];
				C[(i + 1) * r + j + 2] += A[(i + 1) * q + k] * B[k * r + j + 2];
				C[(i + 2) * r + j] += A[(i + 2) * q + k] * B[k * r + j];
				C[(i + 2) * r + j + 1] += A[(i + 2) * q + k] * B[k * r + j + 1];
				C[(i + 2) * r + j + 2] += A[(i + 2) * q + k] * B[k * r + j + 2];
			}
		}
		for(int j = r - rrem; j < r; j++)
		{
			C[i * r + j] = 0.;
			C[(i + 1) * r + j] = 0.;
			C[(i + 2) * r + j] = 0.;
			for(int k = 0; k < q; k++)
			{
				C[i * r + j] += A[i * q + k] * B[k * r + j];
				C[(i + 1) * r + j] += A[(i + 1) * q + k] * B[k * r + j];
				C[(i + 2) * r + j] += A[(i + 2) * q + k] * B[k * r + j];
			}
		}
	}
	for(int i = p - prem; i < p; i++)
	{
		for(int j = 0; j < r - 2; j += 3)
		{
			C[i * r + j] = 0.;
			C[i * r + j + 1] = 0.;
			C[i * r + j + 2] = 0.;
			for(int k = 0; k < q; k++)
			{
				C[i * r + j] += A[i * q + k] * B[k * r + j];
				C[i * r + j + 1] += A[i * q + k] * B[k * r + j + 1];
				C[i * r + j + 2]+= A[i * q + k] * B[k * r + j + 2];
			}
		}
		for(int j = r - rrem; j < r; j++)
		{
			C[i * r + j] = 0.;
			for(int k = 0; k < q; k++)
			{
				C[i * r + j] += A[i * q + k] * B[k * r + j];
			}
		}
	}
	return 0;
}






static inline double ANorm_MPI(double *A, double *F, int N, int m, int p, int t)
{
	int k = N / m, l = N - k * m;
	double* block = nullptr;
	double locMax = -1., glMax = -1.;
	MPI_Comm world = MPI_COMM_WORLD;


	for (int blockj = 0; blockj < procawidthfullonly(N, m, p, t); blockj++)
	{
		for (int j = 0; j < m; j += 3)
		{
			F[j] = 0.;
			F[j + 1] = 0.;
			F[j + 2] = 0.;
		}

		for (int blocki = 0; blocki < k; blocki++)
		{
			block = getBlock(blocki, blockj, N, awidth(N, m, p) * m, m, A);
			// #ifdef __DEBUG__		
			// if (t == 0)
			// {
			// 	printf("NB (%d, %d)\n", blocki, blockj);
			// 	display(block, m, m, m, m);
			// }
			// #endif
			for (int i = 0; i < m; i += 3)
			{
				for (int j = 0; j < m; j +=3)
				{
					F[j] += fabs(block[i * m + j]) + fabs(block[(i + 1) * m + j]) + fabs(block[(i + 2) * m + j]);
					F[j + 1] += fabs(block[i * m + j + 1]) + fabs(block[(i + 1) * m + j + 1]) + fabs(block[(i + 2) * m + j + 1]);
					F[j + 2] += fabs(block[i * m + j + 2]) + fabs(block[(i + 1) * m + j + 2]) + fabs(block[(i + 2) * m + j + 2]);
				}
			}
			// #ifdef __DEBUG__		
			// if (t == 0)
			// {
			// 	printf("F : \n");
			// 	for (int prt = 0; prt < m; prt++)
			// 	{
			// 		printf("%lf  |  ", F[prt]);
			// 	}
			// 	printf("\n");
			// }
			// #endif
		}
		if (l != 0)
		{
			block = getBlock(k, blockj, N, awidth(N, m, p) * m, m, A);
			// #ifdef __DEBUG__		
			// if (t == 0)
			// {
			// 	printf("NB (%d, %d)\n", k, blockj);
			// 	display(block, l, m, m, m);
			// }
			// #endif
			for (int i = 0; i < l; i++)
			{
				for (int j = 0; j < m; j += 3)
				{
					F[j] += fabs(block[i * m + j]);
					F[j + 1] += fabs(block[i * m + j + 1]);
					F[j + 2] += fabs(block[i * m + j + 2]);
				}
			}
			// #ifdef __DEBUG__		
			// if (t == 0)
			// {
			// 	printf("F : \n");
			// 	for (int prt = 0; prt < m; prt++)
			// 	{
			// 		printf("%lf  |  ", F[prt]);
			// 	}
			// 	printf("\n");
			// }
			// #endif
		}

		for (int j = 0; j < m; j += 3)
		{
			locMax = (locMax > F[j] ? locMax : F[j]);
			locMax = (locMax > F[j + 1] ? locMax : F[j + 1]);
			locMax = (locMax > F[j + 2] ? locMax : F[j + 2]);
		}
	}
	if ((t == (p - 1)) && (l != 0))
	{
		for (int j = 0; j < l; j++)
		{
			F[j] = 0.;
		}

		for (int blocki = 0; blocki < k; blocki++)
		{
			block = getBlock(blocki, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);
			// #ifdef __DEBUG__		
			// if (t == 0)
			// {
			// 	printf("NB (%d, %d)\n", blocki, k);
			// 	display(block, m, l, m, m);
			// 	printf("A : \n");
			// 	for (int prt = 0; prt < m * l; prt++)
			// 	{
			// 		printf("%lf  |  ", block[prt]);
			// 	}
			// 	printf("\n");

			// }
			// #endif
			for (int i = 0; i < m; i += 3)
			{
				for (int j = 0; j < l; j++)
				{
					F[j] += fabs(block[i * l + j]) + fabs(block[(i + 1) * l + j]) + fabs(block[(i + 2) * l + j]);
				}
			}
			// #ifdef __DEBUG__		
			// if (t == 0)
			// {
			// 	printf("F : \n");
			// 	for (int prt = 0; prt < l; prt++)
			// 	{
			// 		printf("%lf  |  ", F[prt]);
			// 	}
			// 	printf("\n");
			// }
			// #endif
		}
		block = getBlock(k, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);
		// #ifdef __DEBUG__		
		// if (t == 0)
		// {
		// 	printf("NB (%d, %d)\n", k, k);
		// 	display(block, l, l, l, l);
		// }
		// #endif
		for (int i = 0; i < l; i++)
		{
			for (int j = 0; j < l; j++)
			{
				F[j] += fabs(block[i * l + j]);
			}
		}
		// #ifdef __DEBUG__		
		// if (t == 0)
		// {
		// 	printf("F : \n");
		// 	for (int prt = 0; prt < l; prt++)
		// 	{
		// 		printf("%lf  |  ", F[prt]);
		// 	}
		// 	printf("\n");
		// }
		// #endif

		for (int j = 0; j < l; j++)
		{
			locMax = (locMax > F[j] ? locMax : F[j]);
		}
	}

	MPI_Allreduce(&locMax, &glMax, 1, MPI_DOUBLE, MPI_MAX, world);

	return glMax;
} 


static inline double blockNorm(double *A, int height, int width, int m, double *F)
{
	int kh = height / m, lh = height - kh * m;
	int kw = width / m, lw = width - kw * m;
	int blocki, blockj;
	int ii, jj;
	double *block;
	double ret = -1.;


	//debug
	// printf("AAA %p\n", A);
	// display(A, height, width, m, height);
	//

	for (blockj = 0; blockj < kw; blockj++)
	{
		for (jj = 0; jj < m; jj++)
			F[jj] = 0.;
		for (blocki = 0; blocki < kh; blocki++)
		{
			block = getBlock(blocki, blockj, height, width, m, A);
			for (ii = 0; ii < m; ii++)
			{
				for (jj = 0; jj < m; jj++)
				{
					F[jj] += fabs(block[ii * m + jj]);
				}
			}
		}
		if (lh != 0)
		{
			block = getBlock(kh, blockj, height, width, m, A);
			for (ii = 0; ii < lh; ii++)
			{
				for (jj = 0; jj < m; jj++)
				{
					F[jj] += fabs(block[ii * m + jj]);
				}
			}
		}
		for (jj = 0; jj < m; jj++)
		{
			ret = (ret > F[jj] ? ret : F[jj]);
		}
	}
	if (lw != 0)
	{
		for (jj = 0; jj < lw; jj++)
			F[jj] = 0.;
		for (blocki = 0; blocki < kh; blocki++)
		{
			block = getBlock(blocki, kw, height, width, m, A);
			for (ii = 0; ii < m; ii++)
			{
				for (jj = 0; jj < lw; jj++)
				{
					F[jj] += fabs(block[ii * lw + jj]);
				}
			}
		}
		if (lh != 0)
		{
			block = getBlock(kh, kw, height, width, m, A);
			for (ii = 0; ii < lh; ii++)
			{
				for (jj = 0; jj < lw; jj++)
				{
					F[jj] += fabs(block[ii * lw + jj]);
				}
			}
		}
		for (jj = 0; jj < lw; jj++)
		{
			ret = (ret > F[jj] ? ret : F[jj]);
		}
	}
	return ret;
}


static inline int calculateResidual(double *A, double *B, double *X, double *arg_F, double *dmmt, int N, int m, int p, int t, double scale, double *res)
{
	int k = N / m, l = N - k * m;
	double norm_b = 0., norm_res = 0.;
	double *Ablock = nullptr, *Bblock = nullptr, *Fblock = nullptr, *Xblock = nullptr, *Tempblock = getBlock(0, 0, N, m, m, arg_F), *F = getBlock(1, 0, N, m, m, arg_F);
	MPI_Comm world = MPI_COMM_WORLD;


	MPI_Bcast(X, bsize(N, m), MPI_DOUBLE, printproc, world);

	if( t == printproc)
		norm_b = blockNorm(B, N, 1, m, Tempblock);

	for (int blocki = 0; blocki < k; blocki++)
	{
		Fblock = getBlock(blocki, 0, N, 1, m, F);
		for (int ii = 0; ii < m; ii += 3)
		{
			Fblock[ii] = 0.;
			Fblock[ii + 1] = 0.;
			Fblock[ii + 2] = 0.;
		}
	}
	if (l != 0)
	{
		Fblock = getBlock(k, 0, N, 1, m, F);
		for (int ii = 0; ii < l; ii++)
			Fblock[ii] = 0.;
	}

	for (int blockj = 0; blockj < procawidthfullonly(N, m, p, t); blockj++)
	{
		Xblock = getBlock(L2GBJ(blockj, N, m, p, t), 0, N, 1, m, X);
		for (int blocki = 0; blocki < k; blocki++)
		{
			Ablock = getBlock(blocki, blockj, N, awidth(N, m, p) * m, m, A);
			Fblock = getBlock(blocki, 0, N, 1, m, F);
			matrixMultFrag(Ablock, Xblock, Tempblock, m, m, 1);
			//debug
			// if (t == (p - 1))
			// {
			// printf("-----------------%d :  (%d, %d)-----------------\n", t, blocki, L2GBJ(blockj, N, m, p, t));
			// printf("X block:\n");
			// display(Xblock, m, 1, m, m);
			// printf("A block:\n");
			// display(Ablock, m, m, m, m);
			// printf("Multiplication result:\n");
			// display(Tempblock, m, 1, m, m);
			// printf("Accumulator block\n");
			// display(Fblock, m, 1, m, m);
			// printf("=================%d :  (%d, %d)=================\n", t, blocki, L2GBJ(blockj, N, m, p, t));
			// }
			//
			for (int ii = 0; ii < m; ii += 3)
			{
				Fblock[ii] += Tempblock[ii];
				Fblock[ii + 1] += Tempblock[ii + 1];
				Fblock[ii + 2] += Tempblock[ii + 2];
			}
		}
		if (l != 0)
		{
			Ablock = getBlock(k, blockj, N, awidth(N, m, p) * m, m, A);
			Fblock = getBlock(k, 0, N, 1, m, F);
			matrixMultFrag(Ablock, Xblock, Tempblock, l, m, 1);
			//debug
			// if (t == (p - 1))
			// {
			// printf("-----------------%d :  (%d, %d)-----------------\n", t, k, blockj);
			// printf("X block:\n");
			// display(Xblock, m, 1, m, m);
			// printf("A block:\n");
			// display(Ablock, l, m, m, m);
			// printf("Multiplication result:\n");
			// display(Tempblock, l, 1, l, l);
			// printf("Accumulator block\n");
			// display(Fblock, l, 1, l, l);
			// printf("=================%d :  (%d, %d)=================\n", t, k, blockj);
			// }
			//
			for (int ii = 0; ii < l; ii++)
				Fblock[ii] += Tempblock[ii];
		}
	}



	if ((t == (p - 1)) && (l != 0))
	{
		Xblock = getBlock(k, 0, N, 1, m, X);
		for (int blocki = 0; blocki < k; blocki++)
		{
			Ablock = getBlock(blocki, procawidth(N, m, p, t) - 1, N, awidth(N, m, p) * m, m, A);
			Fblock = getBlock(blocki, 0, N, 1, m, F);
			matrixMultFrag(Ablock, Xblock, Tempblock, m, l, 1);
			//debug
			//if (t == 0)
				// printf("lCheck\n");
			//
			//debug
			// printf("XXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
			// display(Xblock, l, 1, l, l);
			// printf("OOOOOOOOOOO\n");
			// display(Ablock, m, l, m, m);
			// printf("AAAAAAAAAAA\n");
			// display(Tempblock, m, 1, m, m);
			// printf("BBBBBBBBBBB\n");
			// display(Fblock, m, 1, m, m);
			//
			for (int ii = 0; ii < m; ii += 3)
			{
				Fblock[ii] += Tempblock[ii];
				Fblock[ii + 1] += Tempblock[ii + 1];
				Fblock[ii + 2] += Tempblock[ii + 2];
			}
		}
		Ablock = getBlock(k, procawidth(N, m, p, t) - 1, N, awidth(N, m, p) * m, m, A);
		Fblock = getBlock(k, 0, N, 1, m, F);
		matrixMultFrag(Ablock, Xblock, Tempblock, l, l, 1);
		//debug
		// printf("XXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		// display(Xblock, l, 1, l, l);
		// printf("OOOOOOOOOOO\n");
		// display(Ablock, l, l, l, l);
		// printf("AAAAAAAAAAA\n");
		// display(Tempblock, l, 1, l, l);
		// printf("BBBBBBBBBBB\n");
		// display(Fblock, l, 1, l, l);
		//
		for (int ii = 0; ii < l; ii++)
			Fblock[ii] += Tempblock[ii];
	}


	MPI_Allreduce(F, dmmt, bsize(N, m), MPI_DOUBLE, MPI_SUM, world);


	if (t == printproc)
	{
		for (int blocki = 0; blocki < k; blocki++)
		{
			Bblock = getBlock(blocki, 0, N, 1, m, B);
			Fblock = getBlock(blocki, 0, N, 1, m, dmmt);
			for (int ii = 0; ii < m; ii += 3)
			{
				Fblock[ii] -= Bblock[ii];
				Fblock[ii + 1] -= Bblock[ii + 1];
				Fblock[ii + 2] -= Bblock[ii + 2];
			}
		}
		if (l != 0)
		{
			Bblock = getBlock(k, 0, N, 1, m, B);
			Fblock = getBlock(k, 0, N, 1, m, dmmt);
			for (int ii = 0; ii < l; ii++)
			{
				Fblock[ii] -= Bblock[ii];
			}
		}

		norm_res = blockNorm(dmmt, N, 1, m, Tempblock);

		if (fabs(norm_b) >= scale * 1.e-16)
			*res = norm_res / norm_b;
		else
		{
			*res = -1.;
			return -1;
		}
	}
	else
	{
		*res = 0.;
	}


	return 0;
}


static inline double calculateError(double* X, int N, int m)
{
	int k = N / m, l = N - k * m;
	double res = 0.0;
	double *block;

	for (int blocki = 0; blocki < k; blocki++)
	{
		block = getBlock(blocki, 0, N, 1, m, X);
		for (int ii = 0; ii < m; ii += 3)
		{
			res += fabs(block[ii] - (double)((((m * blocki) + ii + 1) & 1)));
			res += fabs(block[ii + 1] - (double)((((m * blocki) + ii) & 1)));
			res += fabs(block[ii + 2] - (double)((((m * blocki) + ii + 1) & 1)));
		}
	}
	if (l != 0)
	{
		block = getBlock(k, 0, N, 1, m, X);
		for (int ii = 0; ii < l; ii++)
		{
			res += fabs(block[ii] - (double)((((m * k) + ii + 1) & 1)));
		}
	}

	return res/(double)((N+1)/2);
}

















int main(int argc, char **argv)
{
	int err = 0;
	double *A = nullptr, *B = nullptr, *F = nullptr, *X = nullptr;
	double *proctime = nullptr;

	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

	MPI_Init(&argc, &argv);
	err = main_do(argc, argv, &A, &B, &F, &X, &proctime);

	if (!(err & 1))
		delete[] A;
	if (!(err & 2))
		delete[] B;
	if (!(err & 4))
		delete[] F;
	if (!(err & 8))
		delete[] X;
	if (!(err & 16))
		delete[] proctime;
		
	MPI_Finalize();

	return 0;
}




int main_do(int argc, char **argv, double **arg_A, double **arg_B, double **arg_F, double **arg_X, double **arg_proctime)
{
	double *A = nullptr, *B = nullptr, *F = nullptr, *X = nullptr;
	double *proctime = nullptr;
	//double *block;
	// double *dmmt = nullptr;
	double cof;
	double error = -1., residual = -1.;
	double timer = 0., fulltime = 0.; 
	int p, t;
	int N, m, r, s;
	// int size;
	//char bar;
	int err = 0, ret = 0;
	char* filename = nullptr;
	MPI_Comm world = MPI_COMM_WORLD;
	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &t);



	if ((argc != 6) && (argc != 5))
	{
		if (t == printproc)
			printf("\nBad argc\n Usage: %s N m r s filename\n", argv[0]);
		err = -t;
	}
	if ((!sscanf(argv[1], "%d", &N)) || (N <= 0))
	{
		if (t == printproc)
			printf("\nBad N\n Usage: %s N m r s filename\n", argv[1]);
		err = -t;
	}
	if ((!sscanf(argv[2], "%d", &m)) || (m <= 0))
	{
		if (t == printproc)
			printf("\nBad m\n Usage: %s N m r s filename\n", argv[2]);
		err = -t;
	}
	else if ((m % 3) != 0)
	{
		if (t == printproc)
			printf("\nm shold be divisible by 3: %s\n Usage: %s N m r s filename\n", argv[2], argv[0]);
		err = -t;
	}
	if ((!sscanf(argv[3], "%d", &r)) || (r <= 0))
	{
		if (t == printproc)
			printf("\nBad r\n Usage: %s N m r s filename\n", argv[3]);
		err = -t;
	}
	if ((!sscanf(argv[4], "%d", &s)) || (s < 0))
	{
		if (t == printproc)
			printf("\nBad s\n Usage: %s N m r s filename\n", argv[4]);
		err = -t;
	}
	if (s == 0)
	{
		if (argc != 6)
		{
			if (t == printproc)
				printf("\nFilename missing\n Usage: %s N m r s filename\n", argv[0]);
			err = -t;
		}
		else
		{
			filename = argv[5];
		}
	}
	MPI_Allreduce(&err, &ret, 1, MPI_INT, MPI_MIN, world);
	if (ret)
	{
		return 0;
	}



	err = 0;
	*(arg_A) = new double[asizempi(N, m, p)];
	A = *(arg_A);
	if (nullptr == A)
	{
		err = err | 1;
	}
	*(arg_B) = new double[bsize(N, m)];
	B = *(arg_B);
	if (nullptr == B)
	{
		err = err | 2;
	}
	*(arg_F) = new double[fsize(N, m) + bsize(N, m)];
	F = *(arg_F);
	if (nullptr == F)
	{
		err = err | 4;
	}
	*(arg_X) = new double[bsize(N, m)];
	X = *(arg_X);
	if (nullptr == X)
	{
		err = err | 8;
	}
	*(arg_proctime) = new double[p];
	proctime = *(arg_proctime);
	if (nullptr == proctime)
	{
		err = err | 16;
	}
	MPI_Allreduce(&err, &ret, 1, MPI_INT, MPI_MIN, world);
	if (ret)
	{
		if (t == printproc)
			printf("memory allocation failed");
		return err;
	}


	err = init(A, B, F, N, m, p, t, s, filename);
	MPI_Allreduce(&err, &ret, 1, MPI_INT, MPI_MIN, world);
	if (ret)
	{
		switch(ret)
		{
		case -1:
			if(t == printproc)
				printf("A is NULL\n");
			break;
		case -2:
			if(t == printproc)
				printf("Can't open file\n");
			break;
		case -3:
			if(t == printproc)
				printf("Bad element encountered\n");
			break;
		case -4:
			if(t == printproc)
				printf("Not enough elements\n");
			break;
		default:
			if(t == printproc)
				printf("Unknown matrix initialization error\n");
		}
		return 0;
	}


	memset(X, 0, bsize(N, m) * sizeof(double));

	err = displayA_MPI(A, F, N, m, p, t, r);
	MPI_Allreduce(&err, &ret, 1, MPI_INT, MPI_MIN, world);
	//debug
	// {
	// 	if (t == 0)
	// 		printf("%d\n\n", ret);
	// }
	if (ret)
	{
		if (t == printproc)
			printf("Ne ponyal?\n");
		return 0;
	}

	if (t == printproc)
	{
		if (t != (p - 1))
			MPI_Recv(B, bsize(N, m), MPI_DOUBLE, p - 1, 0, world, &status);
		err = display(B, N, 1, m, r);
	}
	else if ((t == (p - 1)) && (t != printproc))
		MPI_Send(B, bsize(N, m), MPI_DOUBLE, printproc, 0, world); 
	MPI_Allreduce(&err, &ret, 1, MPI_INT, MPI_MIN, world);
	//debug
	// {
	// 	if (t == 0)
	// 		printf("%d\n\n", ret);
	// }
	if (ret)
	{
		if (t == printproc)
			printf("Ne ponyal?\n");
		return 0;
	}


	// #ifdef __DEBUG__
	// printf("F : %p\n", (void*)F);
	// #endif

	cof = ANorm_MPI(A, F, N, m, p, t);
	// printf("Norm of A = %e\n", cof);
	
	// #ifdef __DEBUG__
	// printf("F : %p\nANorm : %lf\n", (void*)F, cof);
	// #endif


	// 	#ifdef __DEBUG__
	// 	if(t == 0)
	// 	{
	// 		printf("A : %p\nB : %p\nF : %p\nX : %p\n", (void*)A, (void*)B, (void*)F, (void*)X);
	// 	}
	// 	#endif
	timer = MPI_Wtime();
	err = solve_MPI(A, F, B, N, m, p, t, cof);
	timer = MPI_Wtime() - timer;
	//timer /= MPI_Wtick();

	if (err)
	{
		if (t == printproc)
		{
			printf("Algorithm is unapplicable\nFailure on %d step\n", -err);
		}
		return 0;
	}


	if (t == printproc)
	{
		// printf("result:\n B:\n");
		// display(B, N, 1, m, r);
		// printf("X:\n");
		memcpy(X, B, bsize(N, m) * sizeof(double));
		display(X, N, 1, m, r);
	}

	MPI_Gather(&timer, 1, MPI_DOUBLE, proctime, 1, MPI_DOUBLE, printproc, world);

	if (t == printproc)
	{
		fulltime = proctime[0];
		for(int i = 0; i < p; i++)
		{
			if (t == printproc)
				printf("%d process time: %.2f\n", i, proctime[i]);
			if (proctime[i] > fulltime)
				fulltime = proctime[i];
		}
	}

	err = init(A, B, F, N, m, p, t, s, filename);

	err = calculateResidual(A, B, X, F, F + fsize(N, m), N, m, p, t, cof, &residual);

	//#ifndef __FULLRES__
	//if ((N <= 100) || (p > 1))
	//{
	//	err = calculateResidual(A, B, X, F, F + fsize(N, m), N, m, p, t, &residual);
	//}
	//else
	//{
	//	residual = -1.;
	//}
	//#else
	//err = calculateResidual(A, B, X, F, F + fsize(N, m), N, m, p, t, &residual);
	//#endif

	#ifdef __RESO__
	if (p != 1)
		residual = 0.;
	#endif


	if (t == printproc)
	{
		printf("\n\n");
		error = calculateError(X, N, m);
		printf("Error : %e\n", error);
		printf ("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d p = %d\n", argv[0], residual, fulltime, s, N, m, p);
	}



	return 0;
}



