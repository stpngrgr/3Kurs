#include <pthread.h>




unsigned int MSR_offD_size(unsigned int nx, unsigned int ny);
int offDMax(unsigned int /*nx*/, unsigned int /*ny*/);
int offD(unsigned int nx, unsigned int ny, unsigned int l, unsigned int *J);
int offDNum(unsigned int nx, unsigned int ny, unsigned int l);
int MSR_offD_fill(unsigned int nx, unsigned int ny, unsigned int l, unsigned int *J,
                  double &ad, double *A);
int fill_I(unsigned int *I, unsigned int nx, unsigned int ny,
           int p, int t);
int fill_MSR(double *A, unsigned int *I, unsigned int nx, unsigned int ny, 
             int p, int t, pthread_barrier_t *barrier);
int fill_pts(unsigned int nx, unsigned int ny, double x1, double x2, double y1, double y2,
             unsigned int i, unsigned int j, double (*f)(double, double), double fMax, int outlier, double *pts);
double B_row(unsigned int nx, unsigned int ny, double x1, double x2, double y1, double y2,
             unsigned int l, double (*f)(double, double), double fMax, int outlier, int &err);
int fill_B(double *B, unsigned int nx, unsigned int ny, double x1, double x2, double y1, double y2,
           double (*f)(double, double), double fMax, int outlier, int p, int t, pthread_barrier_t *barrier);
double B_row_debug(unsigned int nx, unsigned int ny, double x1, double x2, double y1, double y2,
    unsigned int l, double (*f)(double, double), double fMax, int outlier, int &err);
