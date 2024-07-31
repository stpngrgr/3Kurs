#include <stdio.h>
#include <math.h>
#include "Gauss.h"

int init(double* A, int N, int m, int s, char* filename)
{
	int k = N / m;
	int l = N - k * m;
	int err;
	FILE* F;

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
				if(1 != fscanf( F, "%lf", A + ((i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m))))
				{
					err = (!feof(F) ? -31 : -32);
					fclose(F);
					return err;
				}
			}
		fclose(F);
		break;
	case 1:
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
				A[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)] = N - (i > j ? i : j);
		break;
	case 2:
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
				A[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)] = (i > j ? i : j) + 1;
		break;
	case 3:
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
				A[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)] = fabs(i - j);
		break;
	case 4:
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
				A[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)] = 1.0/(i + j + 1);
		break;
	default:
		return -2;
	}
	for(int i = 0; i < N; i++)
	{
		A[N * N + i] = 0.0;
		for(int j = 0; j < N; j += 2)
		{
			A[N * N + i] += A[(i / m) * N * m + (j / m) * ((i / m) < k ? m : l) * m + (i % m) * ((j / m) < k ? m : l) + (j % m)];
		}
	}
	return 0;
}


