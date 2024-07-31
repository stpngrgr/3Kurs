#include "Gauss.h"
#include <string.h>
#include <stdio.h>

int solve(int N, int m, double* A, double* B, double* X, double* Am)
{
	int k = N / m, l = N - m * k, mei, err;
	double nmin, temp, scale = blockNorm(A, N, m);

	//Down
	for(int t = 0; t < k; t++)
	{
		//display(N, N, m, A, N);
		nmin = 0.0;
		mei = -1;
		for(int i = t; i < k; i++)
		{
			//MSelect
			mcpy(Am + 2 * m * m, A + i * N * m + t * m * m, m * m * sizeof(double));
			if(0 != (err = matrixInverse(Am + 2 * m * m, Am + m * m, m, scale)))
				continue;
			//display(m, m, m, A + i * N * m + t * m * m, m);
			//display(m, m, m, Am + m * m, m);
			temp = norm(Am + m * m, m, m);
			if((mei < 0) || (nmin > temp))
			{
				nmin = temp;
				mei = i;
				mcpy(Am, Am + m * m, m * m * sizeof(double));
				//printf("nnnnnnn%d, %10.3e\n", i, nmin);
			}
		}
		if(mei < 0)
		{
			return -t - 1;
		}
		//display(m, m, m, Am, m);
		//printf("///////%d\n", mei);
		
		//RowSwap
		if(mei != t)
		{
			for(int i = 0; i < k; i++)
			{
				mcpy(Am + m * m, A + t * N * m + i * m * m, m * m * sizeof(double));
				mcpy(A + t * N * m + i * m * m, A + mei * N * m + i * m * m, m * m * sizeof(double));
				mcpy(A + mei * N * m + i * m * m, Am + m * m, m * m * sizeof(double));
			}
			if(l != 0)
			{
				mcpy(Am + m * m, A + t * N * m + k * m * m, m * l * sizeof(double));
				mcpy(A + t * N * m + k * m * m, A + mei * N * m + k * m * m, m * l * sizeof(double));
				mcpy(A + mei * N * m + k * m * m, Am + m * m, m * l * sizeof(double));
			}
			mcpy(Am + m * m, B + t * m, m * sizeof(double));
			mcpy(B + t * m, B + mei * m, m * sizeof(double));
			mcpy(B + mei * m, Am + m * m, m * sizeof(double));
		}

		//RowMult
		for(int i = t + 1; i < k; i++)
		{
			matrixMultFrag(Am, A + t * N * m + i * m * m, Am + m * m, m, m, m);
			mcpy(A + t * N * m + i * m * m, Am + m * m, m * m * sizeof(double));
		}
		if(l != 0)
		{
			matrixMultFrag(Am, A + t * N * m + k * m * m, Am + m * m, m, m, l);
			mcpy(A + t * N * m + k * m * m, Am + m * m, m * l * sizeof(double));
		}
		matrixMultFrag(Am, B + t * m, Am + m * m, m, m, 1);
		mcpy(B + t * m, Am + m * m, m * sizeof(double));
		for(int i = 0; i < m; i++)
			for(int j = 0; j < m; j++)
			{
				A[t * N * m + t * m * m + i * m + j] = (i != j ? 0.0 : 1.0);
			}

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
			for(int io = 0; io < m; io++)
				for(int jo = 0; jo < m; jo++)
					A[i * N * m + t * m * m + io * m + jo] = 0.0;
		}
		if(l != 0)
		{
			for(int j = t + 1; j < k; j++)
			{
				matrixMultSubFrag(A + k * N * m + t * m * l, A + t * N * m + j * m * m, A + k * N * m + j * m * l, l, m, m);
			}
			matrixMultSubFrag(A + k * N * m + t * m * l, A + t * N * m + k * m * m, A + k * N * m + k * m * l, l, m, l);
			matrixMultSubFrag(A + k * N * m + t * m * l, B + t * m, B + k * m, l, m, 1);
			for(int io = 0; io < l; io++)
				for(int jo = 0; jo < m; jo++)
					A[k * N * m + t * m * l + io * m + jo] = 0.0;
		}
	}
	//display(N, N, m, A, N);
	if(l != 0)
	{
		if(0 != matrixInverse(A + k * N * m + k * m * l, Am, l, scale))
		{
			return -k - 1;
		}
		matrixMultFrag(Am, B + k * m, Am + m * m, l, l, 1);
		mcpy(B + k * m, Am + m * m, l * sizeof(double));
		for(int i = 0; i < l; i++)
			for(int j = 0; j < l; j++)
				A[k * N * m + k * m * l + i * l + j] = (i != j ? 0.0 : 1.0);
	}
	//display(N, N, m, A, 10);
	
	//Up
	if(l != 0)
	{
		for(int i = 0; i < k; i++)
		{
			matrixMultSubFrag(A + i * N * m + k * m * m, B + k * m, B + i * m, m, l, 1);
		}
	}
	for(int t = k - 1; t > 0; t--)
	{
		for(int i = 0; i < t; i++)
		{
			matrixMultSubFrag(A + i * N * m + t * m * m, B + t * m, B + i * m, m, m, 1);
		}
	}

	mcpy(X, B, N * sizeof(double));
	return 0;
}
