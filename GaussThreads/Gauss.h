#ifndef GAUSS_H_
#define GAUSS_H_

#include <pthread.h>


class thread_args
{
public:
	int N;
	int m;
	double* A;
	double* B;
	double* X;
	double* F;
	int p;
	int t;
	int s;
	int r;
	char* filename;
	unsigned int ret;
	double res;
	double thread_time;
	double full_time;
	double thread_restime;
	double full_restime;
	pthread_barrier_t* barrier;
};


// int matrixMult(double* A, double* B, double* C, int p, int q, int r);
// int matrixMultFrag(double* A, double* B, double* C, int p, int q, int r);
// int matrixMultRewriteRight(double* A, double* B, double* C, int p, int r);
// int matrixMultSub(double* A, double* B, double* C, int p, int q, int r);
// int matrixMultSubFrag(double* A, double* B, double* C, int p, int q, int r);
// int matrixAdd(double* A, double* B, int p, int q, int sign);
// int matrixInverse(double* A, double* Ai, int p, double scale);
// double calculateResidual(double* A, double* B, double* X, double* F, int N, int m);
// double blockNorm(double* M, int n, int m);
// double norm(double* M, int N, int m, double* F);

//int display(int l, int N, int m, double* A, int r);

int setE(double* M, int m);
int setO(double* M, int n, int m);

void* thread_func(void* vargs);
unsigned int init_thread(int p, int t, int N, int m, double* A, double* B, double* F, int s, char* filename, pthread_barrier_t* barrier);

//unsigned int solve_threads(int N, int m, int t, int p, double mtrNorm, double* A, double* B, double* X, double* F, pthread_barrier_t* barrier);
// double* getBlock(int i, int j, int N, int m, double* A);
//int args_init(thread_args* args, int p, int N, int m, double* A, double* B, double* X, double *F, int s, int r, char* filename, pthread_barrier_t* barrier);
//double* getBlock(int i, int j, int l, int n, int m, double* A);
//int blockSizeD(int rows, int cols);
//int blockSizeB(int rows, int cols);

//void reduce_err(int* ret, pthread_barrier_t* barrier);
//void reduce_erru(unsigned int* ret, pthread_barrier_t* barrier);
//void reduce_max(double* s, pthread_barrier_t* barrier);
//void reduce_min_di(double* d, int* i, pthread_barrier_t* barrier);


//int findCardinal_threads(int p, int t, double* A, double* B, int N, int m, double norm, int s, double* F, pthread_barrier_t* barrier);
//int multRow_threads(int p, int t, double* A, double* B, int N, int m, int s, double* F, pthread_barrier_t* barrier);
//int downStep_threads(double* A, double* B, double* F, int N, int m, int p, int t, int s);
//int lastDownFirstUp_threads(double* A, double* B, int N, int m, int p, int t, double norm, double* F, pthread_barrier_t* barrier);
//int upStep_threads(double* A, double* B, double* F, int N, int m, int p, int t, int s, pthread_barrier_t* barrier);



//double thread_time();
//double full_time();
//int solve(int N, int m, double* A, double* B, double* X, double* Am);
//int init(double* A, int N, int m, int s, char* filename);

static inline int blockSizeD(int rows, int cols)
{
	return rows * cols + ((8 - ((rows * cols) & 7)) & 7);
}

static inline double* getBlock(int i, int j, int l, int n, int m, double* A)
{
	int kw = n / m, lw = n - kw * m;
	int kh = l / m, lh = l - kh * m;
	int fullsize = blockSizeD(m, m), fragsizew = blockSizeD(m, lw), fragsizeh = blockSizeD(m, lh);


	return (A + i * (fullsize * kw + fragsizew) + j * (i < kh ? fullsize : fragsizeh));
}


#endif /* GAUSS_H_ */
