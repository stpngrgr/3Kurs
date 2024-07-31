#include <stdio.h>
#include <math.h>
#include "Eigenvalues.h"

int init(double* A, int N, double* res1, double* res2, double* cof, int s, char* filename)
{
	int err = 0;
	FILE* F;
	double ret1 = 0., ret2 = 0., cot = 0., temp;

	if(A == 0)
		return -1;
	switch(s)
	{
	case 0:
		if((F = fopen(filename, "r")) == nullptr)
		{
			return -30;
		}
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				if(1 != fscanf( F, "%lf", &temp))
				{
					err = (!feof(F) ? -31 : -32);
					fclose(F);
					return err;
				}
				A[i * N + j] = temp;
				if (i == j)
					ret1 += A[i * N + j];
				ret2 += A[i * N + j] * A[i * N + j];
			}
		cot = norm(A, N, N);
		for (int i = 1; i < N; i++)
			for (int j = 0; j < i; j++)
			{
				if (fabs(A[i * N + j] - A[j * N + i]) > e_mash * cot)
					err = -2;
			}
		fclose(F);
		break;
	case 1:
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				A[i * N + j] = N - (i > j ? i : j);
				if(i == j)
					ret1 += A[i * N + j];
				ret2 += A[i * N + j] * A[i * N + j];
			}
		cot = norm(A, N, N);
		break;
	case 2:
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				if (i == j)
				{
					A[i * N + j] = 2;
					ret1 += A[i * N + j];
				}
				else if (abs(i - j) == 1)
					A[i * N + j] = -1;
				else
					A[i * N + j] = 0;
				ret2 += A[i * N + j] * A[i * N + j];
			}
		cot = norm(A, N, N);
		break;
	case 3:
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				if ((i == j) && (i < N - 1))
				{
					A[i * N + j] = 1;
					ret1 += A[i * N + j];
				}
				else if (i == N - 1)
					A[i * N + j] = j + 1;
				else if(j == N - 1)
					A[i * N + j] = i + 1;
				else
					A[i * N + j] = 0;
				ret2 += A[i * N + j] * A[i * N + j];
			}
		ret1 += A[N * N - 1];
		cot = norm(A, N, N);
		break;
	case 4:
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				A[i * N + j] = 1.0/(i + j + 1);
				if(i == j)
					ret1 += A[i * N + j];
				ret2 += A[i * N + j] * A[i * N + j];
			}
		cot = norm(A, N, N);
		break;
	default:
		return -100;
	}
	*res1 = ret1;
	*res2 = sqrt(ret2);
	*cof = cot;
	return err;
}


