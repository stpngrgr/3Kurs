#ifndef GAUSS_H_
#define GAUSS_H_


// int matrixMult(double* A, double* B, double* C, int p, int q, int r);
// int matrixMultFrag(double* A, double* B, double* C, int p, int q, int r);
// int matrixMultRewriteRight(double* A, double* B, double* C, int p, int r);
// int matrixMultSub(double* A, double* B, double* C, int p, int q, int r);
// int matrixMultSubFrag(double* A, double* B, double* C, int p, int q, int r);
// int matrixAdd(double* A, double* B, int p, int q, int sign);
// int matrixInverse(double* A, double* Ai, int p, double scale);
// double norm(double* M, int n, int m);
// double blockNorm(double* M, int N, int m, double* AM);

int display(int l, int N, int m, double* A, int r);

// void mcpy(void* dest, void* from, int size);
// void dcpy(double* dest, double* from, int size);


int init(double* A, int N, int m, int s, char* filename);
int solve(int N, int m, double* A, double* B, double* X, double* Am, double mnorm);


static inline void dcpy(double* dest, double* from, int size)
{
	for(int i = 0; i < size; i++)
	{
		dest[i] = from[i];
	}
}

#endif /* GAUSS_H_ */
