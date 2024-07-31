#include "Gauss.h"
#include <string.h>
#include <stdio.h>
#include <math.h>










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


static inline double norm(double* M, int n, int m)
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







/*
void cardCmp(void *a, void *b, int *len, MPI_Datatype *t)
{
	double *A = (double*) a;
	double *B = (double*) a;
	int m = *((int*)(a));

	for (int i = 0; i < *len; i++)
	{
		if (B[2] > A[2])
		{
			memcpy(b, a, blockSizeB(m, m) + 3);
		}
	}
}
*/



























static inline int downStep_MPI(double *A, double *F, double *B, int N, int m, int p, int t, int s, double cof, MPI_Op *op, MPI_Datatype *type)
{
	int k = N / m, l = N - m * k;
	int err;
	int locCard = -1;
	double locCardNorm = -1., tempNorm;
	double *colBlock = nullptr, *colCpy = nullptr, *resBlock = nullptr, *invBlock = nullptr, *tempBlock = nullptr, *col = getBlock(3, 0, N + 3 * m, m, m, F);
	MPI_Comm world = MPI_COMM_WORLD;


	//temporary
	(void)(op);
	(void)(type);
	struct {double a; int b;} loc, gl;
	//temporaryend

	if (getColumnOwner(N, m, p, s) == t)
	{
		colBlock = getBlock(0, getColumnRank(N, m, p, s), N, awidth(N, m, p) * m, m, A);
		memcpy(col, colBlock, colsize(N, m) * sizeof(double));
			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("nenavist");
			// 	printf("F:\n");
			// 	display(F, N + 3 * m, m, m, N);
			// 	printf("A:\n");
			// 	display(A, N, awidth(N, m, p) * m, m, N);
			// }
			// #endif
	}

	
		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf("AHA %d\n", getColumnOwner(N, m, p, s));
		// 	//display(F, N, m, m, N);
		// }
		// #endif

	MPI_Bcast(col, colsize(N, m), MPI_DOUBLE, getColumnOwner(N, m, p, s), world);

	resBlock = getBlock(0, 0, N + 3 * m, m, m, F);
	invBlock = getBlock(1, 0, N + 3 * m, m, m, F);
	tempBlock = getBlock(2, 0, N + 3 * m, m, m, F);

	// #ifdef __DEBUG__		
	// if (t == printproc)
	// {
	// 	printf("F : %p ---- %p\n", F, F + fsize(N, m) + bsize(N, m));
	// 	printf("resBlock : %p ---- %p\n", resBlock, resBlock + blockSizeD(m, m));
	// 	printf("invBlock : %p ---- %p\n", invBlock, invBlock + blockSizeD(m, m));
	// 	printf("tempBlock : %p ---- %p\n", tempBlock, tempBlock + blockSizeD(m, m));
	// }
	// #endif

	for (int q = s + t; q < k; q += p)
	{

		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf("\n00000000000000AHHFFFFFFFFFFFFFFFFFFFFHA ::\n");
		// }
		// #endif

		colBlock = getBlock(q, 0, N, m, m, col);
		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf("AHA %d\n", q);
		// 	display(colBlock, m, m, m, m);
		// 	display(resBlock, m, m, m, m);
		// }
		// #endif
		memcpy(invBlock, colBlock, blockSizeB(m, m));
		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	display(invBlock, m, m, m, m);
		// }
		// #endif
		err = matrixInverse(invBlock, tempBlock, m, cof);
		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf ("%d\n", err);
		// 	display(tempBlock, m, m, m, m);
		// }
		// #endif
		if (!err)
		{
			tempNorm = norm(tempBlock, m, m);
			if ((locCard < 0) || (tempNorm < locCardNorm))
			{
				locCard = q;
				locCardNorm = tempNorm;
				memcpy(resBlock, tempBlock, blockSizeB(m, m));
			}
		}
	}



	// #ifdef __DEBUG__		
	// if (t == printproc)
	// {
	// 	if (((s / p) + ((s % p) > t ? 1 : 0)) * p + t < k)
	// 	{
	// 		printf("colBlock:\n");
	// 		display(colBlock, m, m, m, m);
	// 	}
	// 	printf("invBlock:\n");
	// 	display(invBlock, m, m, m, m);
	// 	printf("tempBlock:\n");
	// 	display(tempBlock, m, m, m, m);
	// 	printf("resBlock: %d\n", locCard);
	// 	display(resBlock, m, m, m, m);
	// 	printf("%d %d\n", locCard, s);
	// }
	// #endif
	
	// if (locCard >= 0)
	// {
	// 	*((int*)(resBlock - 3)) = m;
	// 	*((int*)(resBlock - 2)) = locCard;
	// 	resBlock[-1] = 1./locCardNorm;
	// }
	// else
	// {
	// 	*((int*)(resBlock - 3)) = m;
	// 	*((int*)(resBlock - 2)) = -1;
	// 	resBlock[-1] = 0.;
	// }

	// MPI_Allreduce(resBlock - 3, invBlock, 1, *type, *op, world);
	// locCard = *((int*)(invBlock + 1));
	// invBlock += 3;


	//temporary
	loc.b = locCard;
	loc.a = (locCard >= 0 ? 1./locCardNorm : 0.);
	MPI_Allreduce(&loc, &gl, 1, MPI_DOUBLE_INT, MPI_MAXLOC, world);
	locCard = gl.b;
	//temporaryend

	// #ifdef __DEBUG__
	// if (t == printproc)
	// 	printf("Card : %d\n", locCard);
	// #endif

	if (locCard == -1)
	{
		return -1;
	}
	

	//temporary
	err = matrixInverse(getBlock(locCard, 0, N, m, m, col), invBlock, m, cof);
	//temporaryend


	// #ifdef __DEBUG__		
	// if (t == printproc)
	// {
	// 	printf("invBlock: %d\n", locCard);
	// 	display(invBlock, m, m, m, m);
	// 	printf("F:\n");
	// 	display(F, N, m, m, N);
	// }
	// #endif

	if (locCard != s)
	{


		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf("%d %d\n", locCard, s);
		// 	printf("FFFFFF %d %d\n", s, getColumnRank(N, m, p, k));
		// 	display(getBlock(s, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A), m, l, m, m);
		// }
		// #endif

		colBlock = getBlock(s, 0, N, m, m, col);
		colCpy = getBlock(locCard, 0, N, m, m, col);
		memcpy(colCpy, colBlock, blockSizeB(m, m));
		
		if (t == (p - 1))
		{
			if (l != 0)
			{
				colBlock = getBlock(s, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);
				colCpy = getBlock(locCard, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);
				memcpy(tempBlock, colBlock, blockSizeB(m, l));
				memcpy(colBlock, colCpy, blockSizeB(m, l));
				memcpy(colCpy, tempBlock, blockSizeB(m, l));
			}

			colBlock = getBlock(s, 0, N, 1, m, B);
			colCpy = getBlock(locCard, 0, N, 1, m, B);
			memcpy(tempBlock, colBlock, blockSizeB(m, 1));
			memcpy(colBlock, colCpy, blockSizeB(m, 1));
			memcpy(colCpy, tempBlock, blockSizeB(m, 1));
		}
	}


	for (int q = (s / p) + ((s % p) >= t ? 1 : 0); q < procawidthfullonly(N, m, p, t); q++)
	{

		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf("\nAHHFFFFFFFFFFFFFFFFFFFFHA ::\n");
		// }
		// #endif

		colBlock = getBlock(s, q, N, awidth(N, m, p) * m, m, A);
		if (locCard != s)
		{
			colCpy = getBlock(locCard, q, N, awidth(N, m, p) * m, m, A);
			memcpy(tempBlock, colBlock, blockSizeB(m, m));
			memcpy(colBlock, colCpy, blockSizeB(m, m));
			memcpy(colCpy, tempBlock, blockSizeB(m, m));
		}
		
		matrixMultRewriteRight(invBlock, colBlock, tempBlock, m, m);

		for (int ii = s + 1; ii < k; ii++)
		{
			resBlock = getBlock(ii, 0, N, m, m, col);
			colCpy = getBlock(ii, q, N, awidth(N, m, p) * m, m, A);
			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("F:\n");
			// 	display(F, N + 3 * m, m, m, N);
			// 	printf("A:\n");
			// 	display(A, N, awidth(N, m, p) * m, m, N);
			// 	printf("resBlock:\n");
			// 	display(resBlock, m, m, m, m);
			// 	printf("colBlock: %d, %d\n", ii, awidth(N, m, p));
			// 	display(colBlock, m, m, m, m);
			// 	printf("colCpy: %d\n", ii);
			// 	display(colCpy, m, m, m, m);
			// }
			// #endif
			matrixMultSubFrag(resBlock, colBlock, colCpy, m, m, m);
			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("colCpy*: %d\n", ii);
			// 	display(colCpy, m, m, m, m);
			// }
			// #endif
		}
		if (l != 0)
		{
			resBlock = getBlock(k, 0, N, m, m, col);
			colCpy = getBlock(k, q, N, awidth(N, m, p) * m, m, A);
			matrixMultSubFrag(resBlock, colBlock, colCpy, l, m, m);
		}
	}

	if (t == (p - 1))
	{
		colBlock = getBlock(s, 0, N, 1, m, B);
		matrixMultRewriteRight(invBlock, colBlock, tempBlock, m, 1);

		for (int ii = s + 1; ii < k; ii++)
		{
			colCpy = getBlock(ii, 0, N, 1, m, B);
			resBlock = getBlock(ii, 0, N, m, m, col);
			matrixMultSubFrag(resBlock, colBlock, colCpy, m, m, 1);
		}
		if (l != 0)
		{
			colCpy = getBlock(k, 0, N, 1, m, B);
			resBlock = getBlock(k, 0, N, m, m, col);
			matrixMultSubFrag(resBlock, colBlock, colCpy, l, m, 1);

			colBlock = getBlock(s, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);

			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("FFFFFF %d %d\n", s, getColumnRank(N, m, p, k));
			// 	display(getBlock(s, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A), m, l, m, m);
			// }
			// #endif

			matrixMultRewriteRight(invBlock, colBlock, tempBlock, m, l);

			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("FFFFFF\n");
			// 	display(colBlock, m, l, m, m);
			// }
			// #endif

			for (int ii = s + 1; ii < k; ii++)
			{
				colCpy = getBlock(ii,  getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);
				resBlock = getBlock(ii, 0, N, m, m, col);


				// #ifdef __DEBUG__		
				// if (t == 0)
				// {
				// 	display(colBlock, m, l, m, m);
				// }
				// #endif


				// #ifdef __DEBUG__		
				// if (t == 0)
				// {
				// 	display(colCpy, m, l, m, m);
				// }
				// #endif


				matrixMultSubFrag(resBlock, colBlock, colCpy, m, m, l);
			}
			colCpy = getBlock(k,  getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);
			resBlock = getBlock(k, 0, N, m, m, col);

			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("Step %d :: \n", s);
			// 	display(colCpy, l, l, l, l);
			// }
			// #endif

			matrixMultSubFrag(resBlock, colBlock, colCpy, l, m, l);

			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("Step %d over :: \n", s);
			// 	display(colCpy, l, l, l, l);
			// }
			// #endif
		}
	}


	// #ifdef __DEBUG__	
	// if (t == printproc)
	// {
	// 	printf("inv :\n");
	// 	display(invBlock, m, m, m, m);
	// 	printf("B :\n");
	// 	display(B, N, 1, m, N);
	// }
	// #endif

	return 0;
}




static inline int lastDownFirstUp_MPI(double *A, double *F, double *B, int N, int m, int p, int t, double cof)
{
	int k = N / m, l = N - k * m;
	double *resBlock = nullptr, *colBlock = nullptr, *tempBlock = nullptr, *col = getBlock(2, 0, N + 3 * m, m, m, F);
	int err = 0;
	MPI_Comm world = MPI_COMM_WORLD;


	if (l != 0)
	{
		if (t == (p - 1))
		{

			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("\n%p ::\n", A);
			// }
			// #endif

			colBlock = getBlock(k, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);

			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("\n%p ::\n", colBlock);
			// }
			// #endif

			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("Last :\n");
			// 	display(colBlock, l, l, l, l);
			// }
			// #endif

			resBlock = getBlock(0, 0, N + 3 * m, m, m, F);
			tempBlock = getBlock(1, 0, N + 3 * m, m, m, F);
			err = matrixInverse(colBlock, resBlock, l, cof);


			// #ifdef __DEBUG__		
			// if (t == printproc)
			// {
			// 	printf("LastInv :\n");
			// 	display(resBlock, l, l, l, l);
			// }
			// #endif

			if (!err)
			{
				
				// #ifdef __DEBUG__	
				// printf("Col :\n");
				// display(F, N, m, m, N);
				// printf("B :\n");
				// display(B, N, 1, m, N);
				// #endif

				colBlock = getBlock(k, 0, N, 1, m, B);
				matrixMultRewriteRight(resBlock, colBlock, tempBlock, l, 1);
				colBlock = getBlock(0, getColumnRank(N, m, p, k), N, awidth(N, m, p) * m, m, A);
				memcpy(col, colBlock, colsize(N, m) * sizeof(double));
				memcpy(col + colsize(N, m) + 1, B, bsize(N, m) * sizeof(double));
				col[colsize(N, m)] = 1.;

				// #ifdef __DEBUG__		
				// if (t == printproc)
				// {
				// 	printf("\nAHHHHHA ::\n");
				// }
				// #endif
			}
			else
			{
				col[colsize(N, m)] = -1.;
			}

			// #ifdef __DEBUG__	
			// printf("Col :\n");
			// display(F, N, m, m, N);
			// printf("B :\n");
			// display(B, N, 1, m, N);
			// #endif



			// #ifdef __DEBUG__	
			// printf("F : \n\n\n");
			// for (int prt = 0; prt < colsize(N, m) + bsize(N, m) + 1; prt++)
			// {
			// 	printf("%lf  |  \n", F[prt]);
			// }
			// printf("\n\n\n");
			// #endif
		}

		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf("\nAHHHHHHHHHA ::\n");
		// }
		// #endif

		MPI_Bcast(col, colsize(N, m) + bsize(N, m) + 1, MPI_DOUBLE, p - 1, world);

		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf("Col 2 :\n");
		// 	display(F, N, m, m, N);
		// }
		// #endif

		if (col[colsize(N, m)] < 0.)
		{
			return -1;
		}
		
		if (t != (p - 1))
		{
			memcpy(B, col + colsize(N, m) + 1, bsize(N, m) * sizeof(double));
		}

		// #ifdef __DEBUG__
		// if (t == 0)
		// {
		// 	printf("B :\n");
		// 	display(B, N, 1, m, N);
		// }
		// #endif

		// #ifdef __DEBUG__		
		// if (t == printproc)
		// {
		// 	printf("\n%p ::\n", A);
		// }
		// #endif

		for (int q = t; q < k; q += p)
		{
			colBlock = getBlock(q, 0, N, m, m, col);
			resBlock = getBlock(q, 0, N, 1, m, B);
			tempBlock = getBlock(k, 0, N, 1, m, B);
			matrixMultSubFrag(colBlock, tempBlock, resBlock, m, l, 1);
		}

		// #ifdef __DEBUG__
		// if (t == 0)
		// {
		// 	printf("B :\n");
		// 	display(B, N, 1, m, N);
		// }
		// #endif
	}
	else
	{
		MPI_Bcast(B, bsize(N, m), MPI_DOUBLE, p - 1, world);
		// #ifdef __DEBUG__
		// if (t == printproc)
		// {
		// 	printf("B :\n");
		// 	display(B, N, 1, m, N);
		// }
		// #endif
	}

	return 0;
}





static inline int upStep_MPI(double *A, double *F, double *B, int N, int m, int p, int t, int s)
{
	//int k = N / m, l = N - m * k;
	double *resBlock = nullptr, *colBlock = nullptr, *tempBlock = nullptr;
	MPI_Comm world = MPI_COMM_WORLD;




	if (t == getColumnOwner(N, m, p, s))
	{
		resBlock = getBlock(0, getColumnRank(N, m, p, s), N, awidth(N, m, p) * m, m, A);
		colBlock = getBlock(s, 0, N, 1, m, B);
		memcpy(F, resBlock, blockSizeB(m, m) * s);
		memcpy(F + blockSizeD(m, m) * s, colBlock, blockSizeB(m, 1));
	}

	MPI_Bcast(F, blockSizeD(m, m) * s + blockSizeD(m, 1), MPI_DOUBLE, getColumnOwner(N, m, p, s), world);


	if ((t == 0) && (t != getColumnOwner(N, m, p, s)))
	{
		colBlock = getBlock(s, 0, N, 1, m, B);
		memcpy(colBlock, F + blockSizeD(m, m) * s, blockSizeB(m, 1));
	}

	// #ifdef __DEBUG__
	// if (t == printproc)
	// {
	// 	printf("B :\n");
	// 	display(B, N, 1, m, N);
	// }
	// #endif

	for (int q = t; q < s; q += p)
	{
		colBlock = getBlock(q, 0, m * s, m, m, F);
		resBlock = getBlock(q, 0, N, 1, m, B);
		tempBlock = F + blockSizeD(m, m) * s;
		matrixMultSubFrag(colBlock, tempBlock, resBlock, m, m, 1);
	}

	return 0;
}




























int solve_MPI(double *A, double *F, double *B, int N, int m, int p, int t, double cof)
{
	int k = N / m;
	int err;
	MPI_Op op;
	MPI_Datatype type;


	// MPI_Op_create(cardCmp, 1, &op);
	// MPI_Type_contiguous(blockSizeD(m, m) + 3, MPI_INT, &type);
	// MPI_Type_commit(&type);


	for (int s = 0; s < k; s++)
	{
		err = downStep_MPI(A, F, B, N, m, p, t, s, cof, &op, &type);
		if (err)
		{
			// MPI_Op_free(&op);
			// MPI_Type_free(&type);
			return -s - 1;
		}

	}

	// #ifdef __DEBUG__
	// if (t == printproc)
	// {
	// 	printf("B almost :\n");
	// 	display(B, N, 1, m, N);
	// }
	// #endif

	err = lastDownFirstUp_MPI(A, F, B, N, m, p, t, cof);

	if (err)
	{
		// MPI_Op_free(&op);
		// MPI_Type_free(&type);
		return -k - 1;
	}


	for(int s = k - 1; s > 0; s--)
	{
		err = upStep_MPI(A, F, B, N, m, p, t, s);

		if (err)
		{
			// MPI_Op_free(&op);
			// MPI_Type_free(&type);
			return s + 100;
		}

	}

		// if (t == 0)
		// 	printf("exited : %d\n", err);
	// MPI_Op_free(&op);
	// MPI_Type_free(&type);
	return 0;
}







