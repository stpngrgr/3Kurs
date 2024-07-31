#include "Gauss.h"
#include <string.h>
#include <stdio.h>
#include <math.h>






//non3mult
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
		for(int j = 0; j < r - rrem; j += 3)
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

static inline int matrixAdd(double* A, double* B, int p, int q, int sign)
{
	if(sign < 0)
	{
		for(int i = 0; i < p; i++)
			for(int j = 0; j < q; j++)
				A[i * q + j] -= B[i * q + j];
	}
	else if(sign > 0)
	{
		for(int i = 0; i < p; i++)
			for(int j = 0; j < q; j++)
				A[i * q + j] += B[i * q + j];
	}
	else
	{
		for(int i = 0; i < p * q; i++)
			A[i] = -1.0 * A[i];
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







int solve(int N, int m, double* A, double* B, double* X, double* Am, double scale)
{
	int k = N / m, l = N - m * k, mei, err;
	double nmin, temp; //scale = blockNorm(A, N, m, Am);
	//double* block;

	//Down
	for(int t = 0; t < k; t++)
	{

		nmin = 0.0;
		mei = -1;
		for(int i = t; i < k; i++)
		{
			//MSelect
			dcpy(Am + 2 * m * m, A + i * N * m + t * m * m, m * m);
			if(0 != (err = matrixInverse(Am + 2 * m * m, Am + m * m, m, scale)))
				continue;
			temp = norm(Am + m * m, m, m);
			if((mei < 0) || (nmin > temp))
			{
				nmin = temp;
				mei = i;
				dcpy(Am, Am + m * m, m * m);
			}
		}
		if(mei < 0)
		{
			return -t - 1;
		}
		

		//RowSwap
		if(mei != t)
		{
			for(int i = t; i < k; i++)
			{
				dcpy(Am + m * m, A + t * N * m + i * m * m, m * m);
				dcpy(A + t * N * m + i * m * m, A + mei * N * m + i * m * m, m * m);
				dcpy(A + mei * N * m + i * m * m, Am + m * m, m * m);
			}
			if(l != 0)
			{
				dcpy(Am + m * m, A + t * N * m + k * m * m, m * l);
				dcpy(A + t * N * m + k * m * m, A + mei * N * m + k * m * m, m * l);
				dcpy(A + mei * N * m + k * m * m, Am + m * m, m * l);
			}
			dcpy(Am + m * m, B + t * m, m);
			dcpy(B + t * m, B + mei * m, m);
			dcpy(B + mei * m, Am + m * m, m);
		}

		//debug
		#ifdef __DEBUG__
			printf("OOOOOOOOOOOOOOOOOOOOO:\n");
			display(N, N, m, A, N);
			display(N, 1, m, B, N);
			display(m, m, m, Am, m);
		#endif

		//RowMult
		for(int i = t + 1; i < k; i++)
		{
			matrixMultRewriteRight(Am, A + t * N * m + i * m * m, Am + m * m, m, m);
		}
		if(l != 0)
		{
			matrixMultRewriteRight(Am, A + t * N * m + k * m * m, Am + m * m, m, l);
		}
		matrixMultRewriteRight(Am, B + t * m, Am + m * m, m, 1);

		//debug
		#ifdef __DEBUG__
			printf("AAAAAAAAAAAAAAAAAAAAA:\n");
			display(N, N, m, A, N);
			display(N, 1, m, B, N);
		 #endif

		//Substr
		for(int i = t + 1; i < k; i++)
		{
			for(int j = t + 1; j < k; j++)
			{
				matrixMultSubFrag(A + i * N * m + t * m * m, A + t * N * m + j * m * m, A + i * N * m + j * m * m, m, m, m);
			}
			if(l != 0)
			{
				matrixMultSubFrag(A + i * N * m + t * m * m, A + t * N * m + k * m * m, A + i * N * m + k * m * m, m, m, l);
			}
			matrixMultSubFrag(A + i * N * m + t * m * m, B + t * m, B + i * m, m, m, 1);
		}
		if(l != 0)
		{
			for(int j = t + 1; j < k; j++)
			{
				matrixMultSubFrag(A + k * N * m + t * m * l, A + t * N * m + j * m * m, A + k * N * m + j * m * l, l, m, m);
			}
			matrixMultSubFrag(A + k * N * m + t * m * l, A + t * N * m + k * m * m, A + k * N * m + k * m * l, l, m, l);
			matrixMultSubFrag(A + k * N * m + t * m * l, B + t * m, B + k * m, l, m, 1);
		}

		//debug
		#ifdef __DEBUG__
			printf("BBBBBBBBBBBBBBBBBBBBB:\n");
			display(N, N, m, A, N);
			display(N, 1, m, B, N);
		#endif
	}
	if(l != 0)
	{
		if(0 != matrixInverse(A + k * N * m + k * m * l, Am, l, scale))
		{
			return -k - 1;
		}
		matrixMultRewriteRight(Am, B + k * m, Am + m * m, l, 1);
	}
	
	//Up
	if(l != 0)
	{
		for(int i = 0; i < k; i++)
		{
			matrixMultSubFrag(A + i * N * m + k * m * m, B + k * m, B + i * m, m, l, 1);
		}
	}
	//debug
	#ifdef __DEBUG__
		display(N, N, m, A, N);
		display(N, 1, m, B, N);
		printf("////////////////////\n");
	#endif

	for(int t = k - 1; t > 0; t--)
	{
		for(int i = 0; i < t; i++)
		{
			matrixMultSubFrag(A + i * N * m + t * m * m, B + t * m, B + i * m, m, m, 1);
		}

		//debug
		#ifdef __DEBUG__
			display(N, N, m, A, N);
			display(N, 1, m, B, N);
		#endif
	}

	dcpy(X, B, N);
	return 0;
}



	////////////////
	// 	if(mei != t)
	// 	{
	// 		dcpy(Am + m * m, A + t * N * m, m * m);
	// 		dcpy(A + t * N * m, A + mei * N * m, m * m);
	// 		dcpy(A + mei * N * m, Am + m * m, m * m);
	// 	}
	// 	for(int i = t + 1; i < k; i++)
	// 	{
	// 		if(mei != t)
	// 		{
	// 			dcpy(Am + m * m, A + t * N * m + i * m * m, m * m);
	// 			dcpy(A + t * N * m + i * m * m, A + mei * N * m + i * m * m, m * m);
	// 			dcpy(A + mei * N * m + i * m * m, Am + m * m, m * m);
	// 		}
	// 		block = A + t * N * m + i * m * m;
	// 		matrixMultRewriteRight(Am, block, Am + m * m, m, m);
	// 		for(int j = t + 1; j < k; j++)
	// 		{
	// 			matrixMultSubFrag(A + j * N * m + t * m * m, block, A + j * N * m + i * m * m, m, m, m);
	// 		}
	// 		if(l != 0)
	// 		{
	// 			matrixMultSubFrag(A + k * N * m + t * m * l, block, A + k * N * m + i * m * l, l, m, m);
	// 		}
	// 	}
	// 	if(l != 0)
	// 	{
	// 		if(mei != t)
	// 		{
	// 			dcpy(Am + m * m, A + t * N * m + k * m * m, m * l);
	// 			dcpy(A + t * N * m + k * m * m, A + mei * N * m + k * m * m, m * l);
	// 			dcpy(A + mei * N * m + k * m * m, Am + m * m, m * l);
	// 		}
	// 		block = A + t * N * m + k * m * m;
	// 		matrixMultRewriteRight(Am, block, Am + m * m, m, l);
	// 		for(int j = t + 1; j < k; j++)
	// 		{
	// 			matrixMultSubFrag(A + j * N * m + t * m * m, block, A + j * N * m + k * m * m, m, m, l);
	// 		}
	// 		matrixMultSubFrag(A + k * N * m + t * m * l, block, A + k * N * m + k * m * l, l, m, l);
	// 	}
	// 	if(mei != t)
	// 	{
	// 		dcpy(Am + m * m, B + t * m, m);
	// 		dcpy(B + t * m, B + mei * m, m);
	// 		dcpy(B + mei * m, Am + m * m, m);
	// 	}
	// 	block = B + t * m;
	// 	matrixMultRewriteRight(Am, block, Am + m * m, m, 1);
	// 	for(int j = t + 1; j < k; j++)
	// 	{
	// 		matrixMultSubFrag(A + j * N * m + t * m * m, block, B + j * m, m, m, 1);
	// 	}
	// 	if(l != 0)
	// 	{
	// 		matrixMultSubFrag(A + k * N * m + t * m * l, block, B + k * m, l, m, 1);
	// 	}
	// }
	// if(l != 0)
	// {
	// 	if(0 != matrixInverse(A + k * N * m + k * m * l, Am, l, scale))
	// 	{
	// 		return -k - 1;
	// 	}
	// 	matrixMultRewriteRight(Am, B + k * m, Am + m * m, l, 1);
	// }
	////////////////
