#ifndef GAUSS_H_
#define GAUSS_H_


int matrixMult(double* A, double* B, double* C, int p, int q, int r);
int matrixMultFrag(double* A, double* B, double* C, int p, int q, int r);
int matrixMultSub(double* A, double* B, double* C, int p, int q, int r);
int matrixMultSubFrag(double* A, double* B, double* C, int p, int q, int r);
int matrixAdd(double* A, double* B, int p, int q, int sign);
int matrixInverse(double* A, double* Ai, int p, double scale);
int solve(int N, int m, double* A, double* B, double* X, double* Am);
int init(double* A, int N, int m, int s, char* filename);
double norm(double* M, int n, int m);
double blockNorm(double* M, int N, int m);

int display(int l, int N, int m, double* A, int r);



#endif /* GAUSS_H_ */
