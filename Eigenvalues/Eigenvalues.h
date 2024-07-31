#ifndef EIGENVALUES_H_
#define EIGENVALUES_H_
#define e_mash 1e-16


#include <math.h>

struct node 
{
	int n;
	int ding;
};

int display(int l, int n, double* A, int r);

int init(double* A, int N, double* res1, double* res2, double* cof, int s, char* filename);
int SolveRec(double* A, int N, double e, double* Am, node* list, double cof);
int tridiagonalize(double* A, int N, double* x, double* y, double* z, double cof);





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




#endif