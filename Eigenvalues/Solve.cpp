#include "Eigenvalues.h"
#include <stdio.h>
#include <stdlib.h>


int cmp(const void* a, const void* b)
{
	return (*((const double*) b) > *((const double*) a) ? -1 : 1);
}



int matrixMultFrag(double* A, double* B, double* C, int p, int q, int r)
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







//static inline int matrixMultRewriteRight(double* A, double* B, double* C, int p, int r)
//{
//	int prem = p % 3, rrem = r % 3;
//
//	for(int j = 0; j < r - rrem; j += 3)
//	{
//		for(int i = 0; i < p - prem; i += 3)
//		{
//			C[i] = 0.;
//			C[i + 1] = 0.;
//			C[i + 2] = 0.;
//			C[p + i] = 0.;
//			C[p + i + 1] = 0.;
//			C[p + i + 2] = 0.;
//			C[2 * p + i] = 0.;
//			C[2 * p + i + 1] = 0.;
//			C[2 * p + i + 2] = 0.;
//			for(int k = 0; k < p; k++)
//			{
//				C[i] += A[i * p + k] * B[k * r + j];
//				C[p + i] += A[i * p + k] * B[k * r + j + 1];
//				C[2 * p + i] += A[i * p + k] * B[k * r + j + 2];
//				C[i + 1] += A[(i + 1) * p + k] * B[k * r + j];
//				C[p + i + 1] += A[(i + 1) * p + k] * B[k * r + j + 1];
//				C[2 * p + i + 1] += A[(i + 1) * p + k] * B[k * r + j + 2];
//				C[i + 2] += A[(i + 2) * p + k] * B[k * r + j];
//				C[p + i + 2] += A[(i + 2) * p + k] * B[k * r + j + 1];
//				C[2 * p + i + 2] += A[(i + 2) * p + k] * B[k * r + j + 2];
//			}
//		}
//		for(int i = p - prem; i < p; i++)
//		{
//			C[i] = 0.;
//			C[p + i] = 0.;
//			C[2 * p + i] = 0.;
//			for(int k = 0; k < p; k++)
//			{
//				C[i] += A[i * p + k] * B[k * r + j];
//				C[p + i] += A[(i + 1) * p + k] * B[k * r + j];
//				C[2 * p + i] += A[(i + 2) * p + k] * B[k * r + j];
//			}
//		}
//		for(int i = 0; i < p - prem; i += 3)
//		{
//			B[i * r + j] = C[i];
//			B[i * r + j + 1] = C[p + i];
//			B[i * r + j + 2] = C[2 * p + i];
//			B[(i + 1) * r + j] = C[i + 1];
//			B[(i + 1) * r + j + 1] = C[p + i + 1];
//			B[(i + 1) * r + j + 2] = C[2 * p + i + 1];
//			B[(i + 2) * r + j] = C[i + 2];
//			B[(i + 2) * r + j + 1] = C[p + i + 2];
//			B[(i + 2) * r + j + 2] = C[2 * p + i + 2];
//		}
//		for(int i = p - prem; i < p; i++)
//		{
//			B[i * r + j] = C[i];
//			B[(i + 1) * r + j] = C[p + i];
//			B[(i + 2) * r + j] = C[2 * p + i];
//		}
//	}
//	for(int j = r - rrem; j < r; j++)
//	{
//		for(int i = 0; i < p - prem; i += 3)
//		{
//			C[i] = 0.;
//			C[i + 1] = 0.;
//			C[i + 2] = 0.;
//			for(int k = 0; k < p; k++)
//			{
//				C[i] += A[i * p + k] * B[k * r + j];
//				C[i + 1] += A[(i + 1) * p + k] * B[k * r + j];
//				C[i + 2]+= A[(i + 2) * p + k] * B[k * r + j];
//			}
//		}
//		for(int i = p - prem; i < p; i++)
//		{
//			C[i] = 0.;
//			for(int k = 0; k < p; k++)
//			{
//				C[i] += A[i * p + k] * B[k * r + j];
//			}
//		}
//		for(int i = 0; i < p - prem; i += 3)
//		{
//			B[i * r + j] = C[i];
//			B[(i + 1) * r + j] = C[i + 1];
//			B[(i + 2) * r + j] = C[i + 2];
//		}
//		for(int i = p - prem; i < p; i++)
//		{
//			B[i * r + j] = C[i];
//		}
//	}
//	return 0;
//}
//
////nnon3mult
//static inline int matrixMultSubFrag(double* A, double* B, double* C, int p, int q, int r)
//{
//	double temp[9];
//	int prem = p % 3, rrem = r % 3;
//
//	for(int i = 0; i < p - prem; i += 3)
//	{
//		for(int j = 0; j < r - rrem; j += 3)
//		{
//			temp[0] = 0.;
//			temp[1] = 0.;
//			temp[2] = 0.;
//			temp[3] = 0.;
//			temp[4] = 0.;
//			temp[5] = 0.;
//			temp[6] = 0.;
//			temp[7] = 0.;
//			temp[8] = 0.;
//			for(int k = 0; k < q; k++)
//			{
//				temp[0] += A[i * q + k] * B[k * r + j];
//				temp[1] += A[i * q + k] * B[k * r + j + 1];
//				temp[2] += A[i * q + k] * B[k * r + j + 2];
//				temp[3] += A[(i + 1) * q + k] * B[k * r + j];
//				temp[4] += A[(i + 1) * q + k] * B[k * r + j + 1];
//				temp[5] += A[(i + 1) * q + k] * B[k * r + j + 2];
//				temp[6] += A[(i + 2) * q + k] * B[k * r + j];
//				temp[7] += A[(i + 2) * q + k] * B[k * r + j + 1];
//				temp[8] += A[(i + 2) * q + k] * B[k * r + j + 2];
//			}
//			C[i * r + j] -= temp[0];
//			C[i * r + j + 1] -= temp[1];
//			C[i * r + j + 2] -= temp[2];
//			C[(i + 1) * r + j] -= temp[3];
//			C[(i + 1) * r + j + 1] -= temp[4];
//			C[(i + 1) * r + j + 2] -= temp[5];
//			C[(i + 2) * r + j] -= temp[6];
//			C[(i + 2) * r + j + 1] -= temp[7];
//			C[(i + 2) * r + j + 2] -= temp[8];
//		}
//		for(int j = r - rrem; j < r; j++)
//		{
//			temp[0] = 0.;
//			temp[1] = 0.;
//			temp[2] = 0.;
//			for(int k = 0; k < q; k++)
//			{
//				temp[0] += A[i * q + k] * B[k * r + j];
//				temp[1] += A[(i + 1) * q + k] * B[k * r + j];
//				temp[2] += A[(i + 2) * q + k] * B[k * r + j];
//			}
//			C[i * r + j] -= temp[0];
//			C[(i + 1) * r + j] -= temp[1];
//			C[(i + 2) * r + j] -= temp[2];
//		}
//	}
//	for(int i = p - prem; i < p; i++)
//	{
//		for(int j = 0; j < r - rrem; j += 3)
//		{
//			temp[0] = 0.;
//			temp[1] = 0.;
//			temp[2] = 0.;
//			for(int k = 0; k < q; k++)
//			{
//				temp[0] += A[i * q + k] * B[k * r + j];
//				temp[1] += A[i * q + k] * B[k * r + j + 1];
//				temp[2] += A[i * q + k] * B[k * r + j + 2];
//			}
//			C[i * r + j] -= temp[0];
//			C[i * r + j + 1] -= temp[1];
//			C[i * r + j + 2] -= temp[2];
//		}
//		for(int j = r - rrem; j < r; j++)
//		{
//			temp[0] = 0.;
//			for(int k = 0; k < q; k++)
//			{
//				temp[0] += A[i * q + k] * B[k * r + j];
//			}
//			C[i * r + j] -= temp[0];
//		}
//	}
//	return 0;
//}
//
//static inline int matrixInverse(double* A, double* Ai, int p, double scale)
//{
//	int mei = -1;
//	double maxe = 0.0, temp = 0.0;
//	
//	for(int i = 0; i < p; i++)
//		for(int j = 0; j < p; j++)
//			Ai[i * p + j] = (i != j ? 0.0 : 1.0);
//	//Down
//	for(int i = 0; i < p; i++)
//	{
//		//display(p, p, p, A, p);
//		//display(p, p, p, Ai, p);
//		//Mselect
//		mei = -1;
//		maxe = 0.0;
//		for(int k = i; k < p; k++)
//		{
//			temp = A[k * p + i];
//			if(fabs(temp) < (1e-16) * scale)
//				continue;
//			else if ((mei < 0) || (maxe < fabs(temp)))
//			{
//				mei = k;
//				maxe = temp;
//			}	
//		}
//		if(mei >= 0)
//		{
//			//RowSwap
//			if(mei != i)
//			{
//				for(int j = 0; j < p; j++)
//				{
//					temp = A[mei * p + j];
//					maxe = Ai[mei * p + j];
//					A[mei * p + j] = A[i * p + j];
//					Ai[mei * p + j] = Ai[i * p + j];
//					A[i * p + j] = temp;
//					Ai[i * p + j] = maxe;
//				}
//			}
//			//RowMult
//			for(int j = i + 1; j < p; j++)
//			{
//				A[i * p + j] /= A[i * p + i];
//			}
//			for(int j = 0; j < p; j++)
//				Ai[i * p + j] /= A[i * p + i];
////			A[i * p + i] = 1.0;
//			//Substr
//			for(int j = i + 1; j < p; j++)
//			{
//				for(int k = i + 1; k < p; k++)
//				{
//					A[j * p + k] -= A[j * p + i] * A[i * p + k];
//				}
//				for(int k = 0; k < p; k++)
//					Ai[j * p + k] -= A[j * p + i] * Ai[i * p + k];
////				A[j * p + i] = 0.0;
//			}
//		}
//		else
//		{
//			return -i - 1;
//		}
//	}
//	
//	//Up
//	for(int i = p - 1; i > 0; i--)
//		for(int j = 0; j < i; j++)
//			for(int k = 0; k < p; k++)
//				Ai[j * p + k] -= A[j * p + i] * Ai[i * p + k];
//	//display(p, p, p, Ai, 10);
//	return 0;
//}







int tridiagonalize(double* A, int N, double* x, double* y, double* z, double cof)
{
	double norm, temp, sxy;
	
	for (int i = 0; i < N - 1; i++)
	{
		// #ifdef __DEBUG__
		// printf("AAAAAAA\n");
		// #endif
		norm = 0.;
		temp = 0.;
		for (int j = i + 2; j < N; j++)
		{
			temp += (A[j * N + i] * A[j * N + i]);
			x[j] = A[j * N + i];
		}
		norm = sqrt(temp + (A[(i + 1) * N + i] * A[(i + 1) * N + i]));
		x[i + 1] = A[(i + 1) * N + i] - norm;
//		display(N, 1, x, N);
		temp = sqrt(temp + (x[i + 1] * x[i + 1]));
		if (temp > e_mash * cof)
		{
			for (int j = i + 1; j < N; j++)
			{
				x[j] /= temp;
			}
//			display(N, 1, x, N);
			
			for (int j = i + 1; j < N; j++)
			{
				y[j] = 0.;
				for (int k = i + 1; k < N; k++)
				{
					y[j] += A[j * N + k] * x[k];
				}
			}
//			display(N, 1, y, N);
			
			for (int j = i + 1; j < N; j++)
			{
				sxy = 0.;
				for (int k = i + 1; k < N; k++)
				{
					sxy += x[k] * y[k];
				}
				z[j] = 2 * y[j] - 2 * sxy * x[j];
			}
//			display(N, 1, z, N);
			
	
			A[2 * i] = A[i * N + i];
			A[2 * i + 1] = norm;
			for (int j = i + 1; j < N; j++)
			{
				for (int k = i + 1; k < N; k++)
				{
					A[j * N + k] -= z[j] * x[k] + x[j] * z[k];
				}
			}
		}
		else if ((norm > e_mash * cof) || (fabs(A[i * N + i]) > e_mash * cof))
		{
			A[2 * i] = A[i * N + i];
			A[2 * i + 1] = A[(i + 1) * N + i];
		}
		else
		{
			A[2 * i] = 0.;
			A[2 * i + 1] = 0.;
		}
	}
	A[2 * (N - 1)] = A[N * N - 1];
//	display(N, 2, A, N);
	return 0;
}




int tridstupid(double* A, int N, double* U, double* T, double* x)
{
	double norm, temp;
	
	for (int i = 0; i < N - 2; i++)
	{
		norm = 0.;
		temp = 0.;
		for (int j = i + 2; j < N; j++)
		{
			temp += (A[j * N + i] * A[j * N + i]);
			x[j] = A[j * N + i];
		}
		norm = sqrt(temp + (A[(i + 1) * N + i] * A[(i + 1) * N + i]));
		x[i + 1] = A[(i + 1) * N + i] - norm;
//		display(N, 1, x, N);
		temp = sqrt(temp + (x[i + 1] * x[i + 1]));
		for (int ii = 0; ii < N; ii++)
			for (int jj = 0; jj < N; jj++)
			{
				U[ii * N + jj] = (ii == jj ? 1. : 0.);
			}
		display(1, N - i - 1, x + i + 1, N);
		if (temp > e_mash)
		{
			for (int j = i + 1; j < N; j++)
			{
				x[j] /= temp;
			}
			for (int ii = i + 1; ii < N; ii++)
				for (int jj = i + 1; jj < N; jj++)
				{
					U[ii * N + jj] -= 2 * x[ii] * x[jj];
				}
			matrixMultFrag(U, A, T, N, N, N);
			matrixMultFrag(T, U, A, N, N, N);
		}
		//display
		display(N, N, A, N);
	}
	return 0;
}




// int eigenQuant(double* A, int N, double lambda, double cof, double e, int* nil)
// {
// 	double x, y = 1., alpha = 0., gamma = 1., u, a, b;
// 	double max;
// 	int m = 0, nilcount = 0;
// //	display(N, 2, A, N);
// 	if (N > 1)
// 		for (int i = 0; i < N - 1; i++)
// 		{
// 			if (fabs(A[2 * i] - lambda) > alpha)
// 				alpha = fabs(A[2 * i] - lambda);
// 			if (fabs(A[2 * i + 1]) > alpha)
// 				alpha = fabs(A[2 * i + 1]);
// 		}
// 	if (fabs(A[2 * N - 2] - lambda) > alpha)
// 		alpha = fabs(A[2 * N - 2] - lambda);
// 	alpha *= 4.;
// 	if (fabs(alpha) < e_mash * cof)
// 	{
// 		*nil = N;
// 		return 0;
// 	}
// 	x = (A[0] - lambda) / alpha;
// 	/*if (fabs(x) < e * cof)
// 		nilcount++;
// 	else*/ if (x < 0)
// 		m++;
// //	printf("lambda = %e    step 0:   x = %e,   y = %e\n", lambda, x, y);
// 	for (int i = 1; i < N; i++)
// 	{
// 		a = (A[2 * i] - lambda) / alpha;
// 		b = A[(i - 1) * 2 + 1] / alpha;
// 		max = (fabs(x) > fabs(b * (b * y)) ? fabs(x) : fabs(b * (b * y)));
// 		if (fabs(max) < e * cof)
// 		{
// 			*nil = N - i + nilcount;
// 			return m; 
// 		}
// 		gamma = (1 / e_mash) / max;
// 		u = gamma * ((a * x) - ((b * b) * y));
// //		printf("a = %e, b = %e,u = %e\n", a, b, u);
// 		/*if (fabs(u * x) < e * cof)
// 			nilcount++;
// 		else*/ if (u * x < 0)
// 			m++;
// 		y = gamma * x;
// 		x = u;
// //		printf("lambda = %e    step %d:   x = %e,   y = %e\n", lambda, i, x, y);
// 	}
// 	*nil = nilcount;
// 	return m;
// }







// #ifdef __DEBUG__
// int eigenQuant_(double* A, int N, double lambda, double cof, double e, int rght)
// {
// 	double x, y = 1., alpha = 0., gamma = 1., u, cho, a, b;
// 	double max;
// 	int m = 0;
// 	int yas = 0;

// //	display(N, 2, A, N);
// 	if (N > 1)
// 		for (int i = 0; i < N - 1; i++)
// 		{
// 			if (fabs(A[2 * i] - lambda) > alpha)
// 				alpha = fabs(A[2 * i] - lambda);
// 			if (fabs(A[2 * i + 1]) > alpha)
// 				alpha = fabs(A[2 * i + 1]);
// 		}
// 	if (fabs(A[2 * N - 2] - lambda) > alpha)
// 		alpha = fabs(A[2 * N - 2] - lambda);
// 	if (fabs(alpha) < e_mash * cof)
// 	{
// 		return N;
// 	}
// 	alpha *= 4.;
// 	x = (A[0] - lambda) / alpha;
// 	if ((x <= 0.) || (fabs(x) < cof * e))
// 		m++;
// //	printf("lambda = %e    step 0:   x = %e,   y = %e\n", lambda, x, y);
// 	for (int i = 1; i < N; i++)
// 	{
// 		a = (A[2 * i] - lambda) / alpha;
// 		b = A[(i - 1) * 2 + 1] / alpha;
// 		cho = b * b * y;
// 		max = (fabs(x) > fabs(cho) ? fabs(x) : fabs(cho));
// 		if (fabs(max) < e_mash * cof)
// 		{
// 			return -1; 
// 		}
// 		u = ((a * x) - cho) / (e_mash * max);
// 		// if ((u * x <= 0.) || (fabs(x * u) < e * cof))
// 		// 	m++;
// 		if (((x < 0.) && (u > 0.)) || ((x > 0.) && (u < 0.)) || (fabs(u) < e * cof))
// 		 	m++;
// 		y = x / (e_mash * max);
// 		x = u;
// //		printf("lambda = %e    step %d:   x = %e,   y = %e\n", lambda, i, x, y);
// 	}
	
// 	if (m > rght)
// 	{
// 		y = 1.;
// 		gamma = 1.;
// 		printf("m = %d vs right = %d\n for lambda = %.16e\n alpha = %.16e\n", m, rght, lambda, alpha);
// 		x = (A[0] - lambda) / alpha;
// 		if ((x <= 0.) || (fabs(x) < cof * e))
// 			yas = 1;
// 		printf("x = %.16e | yes? = %d\n\n----------------------\n\n", x, yas);
// 		for (int i = 1; i < N; i++)
// 		{
// 			yas = 0;
// 			a = (A[2 * i] - lambda) / alpha;
// 			b = A[(i - 1) * 2 + 1] / alpha;
// 			cho = b * b * y;
// 			max = (fabs(x) > fabs(cho) ? fabs(x) : fabs(cho));
// 			if (fabs(max) < e_mash * cof)
// 			{
// 				return -1; 
// 			}
// 			u = ((a * x) - cho) / (e_mash * max);
// 			printf("a = %.16e \n b = %.16e \n cho = %.16e \n max = %.16e \n u = %.16e \n ", a, b, cho, max, u);
// 			if ((u * x <= 0.) || (fabs(u * x) < e * cof))
// 				yas = 1;
// 			y = x / (e_mash * max);
// 			x = u;
// 			printf("yes? = %d \n x = %.16e \n y = %.16e \n\n----------------------\n\n", yas, x, y);
// 		}
// 	}

// 	return m;
// }
// #endif




int eigenQuant(double* A, int N, double lambda, double cof, double e)
{
	double x, y = 1., alpha = 0., u, cho, a, b;
	double max;
	int m = 0;

//	display(N, 2, A, N);
	if (N > 1)
		for (int i = 0; i < N - 1; i++)
		{
			if (fabs(A[2 * i] - lambda) > alpha)
				alpha = fabs(A[2 * i] - lambda);
			if (fabs(A[2 * i + 1]) > alpha)
				alpha = fabs(A[2 * i + 1]);
		}
	if (fabs(A[2 * N - 2] - lambda) > alpha)
		alpha = fabs(A[2 * N - 2] - lambda);
	if (fabs(alpha) < e_mash * cof)
	{
		return N;
	}
	alpha *= 4.;
	x = (A[0] - lambda) / alpha;
	if ((x <= 0.) || (fabs(x) < e))
		m++;
//	printf("lambda = %e    step 0:   x = %e,   y = %e\n", lambda, x, y);
	for (int i = 1; i < N; i++)
	{
		a = (A[2 * i] - lambda) / alpha;
		b = A[(i - 1) * 2 + 1] / alpha;
		cho = b * b * y;
		max = (fabs(x) > fabs(cho) ? fabs(x) : fabs(cho));
		if (fabs(max) < e_mash)
		{
			printf("%lf\n", alpha);
			return -1; 
		}
		u = ((a * x) - cho) / (e_mash * max);
		// if ((u * x <= 0.) || (fabs(x * u) < e * cof))
		// 	m++;
		if (((x < 0.) && (u > 0.)) || ((x > 0.) && (u < 0.)) || (fabs(u) < e))
		 	m++;
		y = x / (e_mash * max);
		x = u;
//		printf("lambda = %e    step %d:   x = %e,   y = %e\n", lambda, i, x, y);
	}

	return m;
}







// int recStep(double* A, int N, double ai, double bi, double e, int n_ai, int n_bi, double* arr, double cof)
// {
// 	double a = ai, b = bi, c;
// 	int n_a = n_ai, n_b = n_bi, n_c, flg = 0;
	
// 	while ((flg == 0) && (fabs(b - a) > e * cof))
// 	{
// 		c = (a + b) / 2;
// 		n_c = eigenQuant(A, N, c);
// 		if (n_c == n_a)
// 			a = c;
// 		else if (n_c == n_b)
// 			b = c;
// 		else
// 			flg = 1;
// 	}
// 	if (fabs(b - a) > e * cof)
// 	{
// 		if ((n_c - n_a) > 0)
// 			recStep(A, N, a, c, e, n_a, n_c, arr, cof);
// 		if ((n_b - n_c) > 0)
// 			recStep(A, N, c, b, e, n_c, n_b, arr + (n_c - n_a), cof);
// 	}
// 	else
// 	{
// 		c = (a + b) / 2;
// 		for (int i = 0; i < (n_b - n_a); i++)
// 		{
// 			arr[i] = c;
// 		}
// 	}
// 	return n_b - n_a;
// }



int recStep_(double* A, int N, double a, double b, double e, double* arr, node* list, double cof)
{
	double cura = a, curb = b, c;
	//int nil_cura, nil_curb, nil_c;
	int n_cura = eigenQuant(A, N, a, cof, e), n_curb = eigenQuant(A, N, b, cof, e), n_c, count;
	int eigenum = n_curb - n_cura;
	int timer = 0, timer_limit = 100 * eigenum;
	int dings_left = eigenum;
	int flg;
	node *it;
	
	arr[0] = a;
	arr[2 * N - 1] = b;
	list[0] = {eigenum, 0};
	#ifdef __DEBUG__
	for (int i = 1; i < N; i++)
	{
		list[i] = {0, 0};
	}
	#endif
	if ((n_cura < 0) || (n_curb < 0))
		return -2;
	
	while ((dings_left > 0) && (timer < timer_limit))
	{
		// #ifdef __DEBUG__
		// printf("%d, %d\n", dings_left, timer);
		// #endif
		 for (it = list, count = 0; count < eigenum; count += it->n, it += it->n)
		 {
			// #ifdef __DEBUG__
			// for (int i = 0, cnt = 0; (i < N) && (cnt < N) && (minflg == 0); i++)
			// {
			// 	if (list[i].n > 0)
			// 	{
			// 		printf("%d : ||   %10.4e  |  %d  ,  %10.4e  |  %d   || %d\n", i, arr[2 * cnt], eigenQuant_(A, N, arr[2 * cnt]), arr[2 * (cnt + list[i].n) - 1], eigenQuant_(A, N, arr[2 * (cnt + list[i].n) - 1]), list[i].n);
			// 		cnt += list[i].n;
			// 	}
			// 	else if (list[i].n < 0)
			// 	{
			// 		printf("%d : ||   %10.4e  |  %d  ,  %10.4e  |  %d   || %d\n\n", i, arr[2 * cnt], eigenQuant_(A, N, arr[2 * cnt]), arr[2 * cnt + 1], eigenQuant_(A, N, arr[2 * cnt + 1]), list[i].n);
			// 		printf("-----------------------------------------\n\n");
			// 		minflg = 1;
			// 	}
			// }
			// #endif
			 if (it->ding == 0)
			 {
				 n_cura = count;
				 n_curb = count + it->n;
				 cura = arr[2 * n_cura];
				 curb = arr[2 * (n_curb - 1) + 1];
				 if (fabs(curb - cura) > e * cof)
				 {
					 flg = 0;
					 while ((fabs(curb - cura) > e * cof) && (flg == 0))
					 {
						// #ifdef __DEBUG__
						// printf("%lf, %lf\n", cura, curb);
						// #endif
						 c = (curb + cura) / 2.;
						 n_c = eigenQuant(A, N, c, cof, e);
						 if (n_c < 0)
							return -2;
						// #ifdef __DEBUG__
						// printf("%lf, %lf AAAAAAAAAA\n", cura, curb);
						// #endif
						 timer++;
						 if (n_c == n_cura)
							 cura = c;
						 else if (n_c == n_curb)
							 curb = c;
						 else
						 {
							 #ifdef __DEBUG__
							 if ((n_c - n_cura < 0) || (n_curb - n_c < 0))
							 {
								 printf("%.16e, %d | %.16e, %d | %.16e, %d\n", cura, n_cura, c, n_c, curb, n_curb);
							 }
							 #endif
							 *(it) = {n_c - n_cura, 0};
							 *(it + it->n) = {n_curb - n_c, 0};
							 flg = 1;
						 }
					 }
					 arr[2 * n_cura] = cura;
					 arr[2 * (n_curb - 1) + 1] = curb;
					 if (flg != 0)
					 {
						 arr[2 * (n_c - 1) + 1] = c;
						 arr[2 * (n_c)] = c;
						 count += it->n;
						 it += it->n;
					 }
				 }
				 else
				 {
					 it->ding = it->n;
					 dings_left -= it->ding;
				 }
			 }
			// #ifdef __DEBUG__
			// printf("%p AAAAAAAAAAAAAAAAAAAAAAAA\n", (void*)it);
			// #endif
		 }
	}
	
	if(dings_left == 0)
	{
		for(it = list, count = 0; count < eigenum; count += it->ding, it += it->ding)
		{
			cura = (arr[2 * (count + it->n - 1) + 1] + arr[2 * count]) / 2.;
			for(int i = 0; i < it->ding; i++)
			{
				arr[count + i] =  cura;
			}
		}
		return timer;
	}
	else
		return -1;
}




int SolveRec(double* A, int N, double e, double* Am, node* list, double cof)
{
//	int n_a, n_b;
	int ret = 0, res = 0;
	int count = 0;
	double *arr_start = Am;
//	double *U = new double[N * N], *T = new double[N * N], *AA = new double[N * N];
//	for (int i = 0; i < N * N; i++)
//		AA[i] = A[i];
//	if (0 != tridstupid(A, N, U, T, Am))
//		return -1;
//	for (int ii = 0; ii < N - 1; ii++)
//	{
//		U[2 * ii] = A[ii * N + ii];
//		U[2 * ii + 1] = A[ii * N + ii + 1];
//	}
//	U[2 * (N - 1)] = A[(N - 1) * N + N - 1];
//	if (0 != tridiagonalize(A, N, Am, Am + N, Am + 2 * N, cof))
//		return -1;
	//display
//	printf("\n\n");
//	for (int i = 0; i < 2 * N - 1; i += 2)
//		printf(" %10.3e", U[i]);
//	printf("\n");
//	for (int i = 1; i < 2 * N - 1; i += 2)
//		printf(" %10.3e", U[i]);
//	printf("\n\n");
	//
	//display
//	printf("\n\n");
//	for (int i = 0; i < 2 * N - 1; i += 2)
//		printf(" %10.3e", A[i]);
//	printf("\n");
//	for (int i = 1; i < 2 * N - 1; i += 2)
//		printf(" %10.3e", A[i]);
//	printf("\n\n");
	//
//	n_a = eigenQuant(A, N, -cof * 1.1);
//	n_b = eigenQuant(A, N, cof * 1.1);
//	ret = recStep(A, N, -cof * 1.1, cof * 1.1, e, n_a, n_b, Am, cof);
	for (int i = 0; i < N - 1; i++)
	{
		if (fabs(A[2 * i + 1]) < e_mash * cof)
		{
			if (count != 0)
			{
				#ifdef __DEBUG__
				printf("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa ---------- %d, %d\n", 2 * (i - count), count + 1);
				#endif
				res = recStep_(A + 2 * (i - count), count + 1, -cof * 1.01, cof * 1.01, e, arr_start, list, cof);
				arr_start += count + 1;
				if (res < 0)
					return res;
				else
					ret += res;
			}
			else
			{
				*arr_start = A[2 * i];
				arr_start++;
			}
			count = 0;
		}
		else
			count++;
	}
	if (count != 0)
	{
		#ifdef __DEBUG__
		printf("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa ---------- %d, %d\n", 2 * (N - 1 - count), count + 1);
		#endif
		res = recStep_(A + 2 * (N - 1 - count), count + 1, -cof * 1.01, cof * 1.01, e, arr_start, list, cof);
		if (res < 0)
			return res;
		else
			ret += res;
	}
	else
	{
		*arr_start = A[2 * (N - 1)];
	}
	qsort(Am, N, sizeof(double), cmp);
//	delete[] U;
//	delete[] T;
//	delete[] AA;
	return ret;
}

