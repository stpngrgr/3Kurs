#include "Solve.h"
#include <stdio.h>
#include "fill_matrix.h"







static inline void reduce_err(int* ret, pthread_barrier_t* barrier)
{
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static int global = 0;
	int barret;


	//lock
	pthread_mutex_lock(&m);
	if((global == 0) && (*ret != 0))
		global = *ret;
	pthread_mutex_unlock(&m);
	pthread_barrier_wait(barrier);
	*ret = global;
	barret = pthread_barrier_wait(barrier);
	if(barret == PTHREAD_BARRIER_SERIAL_THREAD)
	{
		global = 0;
	}
	pthread_barrier_wait(barrier);
	return;
}






unsigned int MSR_offD_size(unsigned int nx, unsigned int ny)
{
    return 6 * (nx - 1) * (ny - 1) + 8 * (nx + ny - 2) + 10;
}





static inline void ij2l(unsigned int nx, unsigned int /*ny*/, 
          unsigned int i, unsigned int j, unsigned int &l)
{
    l = i + j * (nx + 1);
    return;
}


static inline void l2ij(unsigned int nx, unsigned int /*ny*/, 
          unsigned int &i, unsigned int &j, unsigned int l)
{
    j = l / (nx + 1);
    i = l - j * (nx + 1);
    return;
}


int offDMax(unsigned int /*nx*/, unsigned int /*ny*/)
{
    return 6;
}


int offD(unsigned int nx, unsigned int ny, unsigned int l, unsigned int *J)
{
    unsigned int i, j;
    l2ij(nx, ny, i, j, l);

    if ( ((i > 0) && (i < nx)) && ((j > 0) && (j < ny)) )
    {
        ij2l(nx, ny, i + 1, j, J[0]);
        ij2l(nx, ny, i, j - 1, J[1]);
        ij2l(nx, ny, i - 1, j - 1, J[2]);
        ij2l(nx, ny, i - 1, j, J[3]);
        ij2l(nx, ny, i, j + 1, J[4]);
        ij2l(nx, ny, i + 1, j + 1, J[5]);
        return 6;
    }


    if ( (i == 0) && ((j > 0) && (j < ny)) )
    {
        ij2l(nx, ny, 1, j, J[0]);
        ij2l(nx, ny, 0, j - 1, J[1]);
        ij2l(nx, ny, 0, j + 1, J[2]);
        ij2l(nx, ny, 1, j + 1, J[3]);
        return 4;
    }

    if ( (i == nx) && ((j > 0) && (j < ny)) )
    {
        ij2l(nx, ny, nx, j - 1, J[0]);
        ij2l(nx, ny, nx - 1, j - 1, J[1]);
        ij2l(nx, ny, nx - 1, j, J[2]);
        ij2l(nx, ny, nx, j + 1, J[3]);
        return 4;
    }

    if ( (j == 0) && ((i > 0) && (i < nx)) )
    {
        ij2l(nx, ny, i + 1, 0, J[0]);
        ij2l(nx, ny, i - 1, 0, J[1]);
        ij2l(nx, ny, i, 1, J[2]);
        ij2l(nx, ny, i + 1, 1, J[3]);
        return 4;
    }

    if ( (j == ny) && ((i > 0) && (i < nx)) )
    {
        ij2l(nx, ny, i + 1, ny, J[0]);
        ij2l(nx, ny, i, ny - 1, J[1]);
        ij2l(nx, ny, i - 1, ny - 1, J[2]);
        ij2l(nx, ny, i - 1, ny, J[3]);
        return 4;
    }


    if ( (i == 0) && (j == 0) )
    {
        ij2l(nx, ny, 1, 0, J[0]);
        ij2l(nx, ny, 0, 1, J[1]);
        ij2l(nx, ny, 1, 1, J[2]);
        return 3;
    }

    if ( (i == nx) && (j == ny) )
    {
        ij2l(nx, ny, nx, ny - 1, J[0]);
        ij2l(nx, ny, nx - 1, ny - 1, J[1]);
        ij2l(nx, ny, nx - 1, ny, J[2]);
        return 3;
    }


    if ( (i == 0) && (j == ny) )
    {
        ij2l(nx, ny, 1, ny, J[0]);
        ij2l(nx, ny, 0, ny - 1, J[1]);
        return 2;
    }

    if ( (i == nx) && (j == 0) )
    {
        ij2l(nx, ny, nx - 1, 0, J[0]);
        ij2l(nx, ny, nx, 1, J[1]);
        return 2;
    }


    return -1;
}


int offDNum(unsigned int nx, unsigned int ny, unsigned int l)
{
    unsigned int i, j;
    l2ij(nx, ny, i, j, l);

    if ( ((i > 0) && (i < nx)) && ((j > 0) && (j < ny)) )
    {
        return 6;
    }


    if ( ( ((i == 0) || (i == nx)) && ((j > 0) && (j < ny)) ) 
        || ( ((j == 0) || (j == ny)) && ((i > 0) && (i < nx)) ) )
    {
        return 4;
    }


    if ( ( ((i == 0)) && (j == 0) ) || ( (i == nx) && (j == ny) ) )
    {
        return 3;
    }

    if ( ( (i == 0) && (j == ny) ) || ( (i == nx) && (j == 0) ) )
    {
        return 2;
    }


    return -1;
}


int MSR_offD_fill(unsigned int nx, unsigned int ny, unsigned int l, unsigned int *J,
                  double &ad, double *A)
{
    unsigned int i, j;
    l2ij(nx, ny, i, j, l);

    if ( ((i > 0) && (i < nx)) && ((j > 0) && (j < ny)) )
    {
        ij2l(nx, ny, i + 1, j, J[0]);
        ij2l(nx, ny, i, j - 1, J[1]);
        ij2l(nx, ny, i - 1, j - 1, J[2]);
        ij2l(nx, ny, i - 1, j, J[3]);
        ij2l(nx, ny, i, j + 1, J[4]);
        ij2l(nx, ny, i + 1, j + 1, J[5]);
        ad = 0.5;
        for (int k = 0; k < 6; k++)
        {
            A[k] = 1. / 12.;
        }
        return 6;
    }


    if ( (i == 0) && ((j > 0) && (j < ny)) ) //left
    {
        ij2l(nx, ny, 1, j, J[0]);
        ij2l(nx, ny, 0, j - 1, J[1]);
        ij2l(nx, ny, 0, j + 1, J[2]);
        ij2l(nx, ny, 1, j + 1, J[3]);
        ad = 0.25;
        A[0] = 1. / 12.;
        A[1] = 1. / 24.;
        A[2] = 1. / 24.;
        A[3] = 1. / 12.;
        return 4;
    }

    if ( (i == nx) && ((j > 0) && (j < ny)) ) //right
    {
        ij2l(nx, ny, nx, j - 1, J[0]);
        ij2l(nx, ny, nx - 1, j - 1, J[1]);
        ij2l(nx, ny, nx - 1, j, J[2]);
        ij2l(nx, ny, nx, j + 1, J[3]);
        ad = 0.25;
        A[0] = 1. / 24.;
        A[1] = 1. / 12.;
        A[2] = 1. / 12.;
        A[3] = 1. / 24.;
        return 4;
    }

    if ( (j == 0) && ((i > 0) && (i < nx)) ) //down
    {
        ij2l(nx, ny, i + 1, 0, J[0]);
        ij2l(nx, ny, i - 1, 0, J[1]);
        ij2l(nx, ny, i, 1, J[2]);
        ij2l(nx, ny, i + 1, 1, J[3]);
        ad = 0.25;
        A[0] = 1. / 24.;
        A[1] = 1. / 24.;
        A[2] = 1. / 12.;
        A[3] = 1. / 12.;
        return 4;
    }

    if ( (j == ny) && ((i > 0) && (i < nx)) ) //up
    {
        ij2l(nx, ny, i + 1, ny, J[0]);
        ij2l(nx, ny, i, ny - 1, J[1]);
        ij2l(nx, ny, i - 1, ny - 1, J[2]);
        ij2l(nx, ny, i - 1, ny, J[3]);
        ad = 0.25;
        A[0] = 1. / 24.;
        A[1] = 1. / 12.;
        A[2] = 1. / 12.;
        A[3] = 1. / 24.;
        return 4;
    }


    if ( (i == 0) && (j == 0) )
    {
        ij2l(nx, ny, 1, 0, J[0]);
        ij2l(nx, ny, 0, 1, J[1]);
        ij2l(nx, ny, 1, 1, J[2]);
        ad = 1. / 6.;
        A[0] = 1. / 24.;
        A[1] = 1. / 24.;
        A[2] = 1. / 12.;
        return 3;
    }

    if ( (i == nx) && (j == ny) )
    {
        ij2l(nx, ny, nx, ny - 1, J[0]);
        ij2l(nx, ny, nx - 1, ny - 1, J[1]);
        ij2l(nx, ny, nx - 1, ny, J[2]);
        ad = 1. / 6.;
        A[0] = 1. / 24.;
        A[1] = 1. / 12.;
        A[2] = 1. / 24.;
        return 3;
    }


    if ( (i == 0) && (j == ny) )
    {
        ij2l(nx, ny, 1, ny, J[0]);
        ij2l(nx, ny, 0, ny - 1, J[1]);
        ad = 1. / 12.;
        A[0] = 1. / 24.;
        A[1] = 1. / 24.;
        return 2;
    }

    if ( (i == nx) && (j == 0) )
    {
        ij2l(nx, ny, nx - 1, 0, J[0]);
        ij2l(nx, ny, nx, 1, J[1]);
        ad = 1. / 12.;
        A[0] = 1. / 24.;
        A[1] = 1. / 24.;
        return 2;
    }


    return -1;
}





int fill_I(unsigned int *I, unsigned int nx, unsigned int ny,
           int p, int t, pthread_barrier_t *barrier)
{
    unsigned int size = MSR_offD_size(nx, ny);
    unsigned int n = (nx + 1) * (ny + 1);
    // unsigned int i, j;
    unsigned int l;
    unsigned int s;
    unsigned int *J, *temp;
    unsigned int r = n;
    int err = 0;


    J = new unsigned int[offDMax(nx, ny)];
    if (!J)
        err = -1;
    //reduce err
    reduce_err(&err, barrier);

    if (!err)
    {
        if (t == 0)
        {
            for (l = 0; l < n; l++)
            {
                s = offDNum(nx, ny, l);
                I[l] = r;
                r += s;
            }
            I[l] = r;
            if(r != n + size + 1)
                err = -1;
        }
        //reduce err 
        reduce_err(&err, barrier);
    }

    if (!err)
    {
        for (l = tRowsBeg(n, p, t); l < tRowsEnd(n, p, t); l++)
        {
            s = offD(nx, ny, l, J);
            if (s != I[l + 1] - I[l])
                err = -1;
            temp = I + I[l];
            for (unsigned int r = 0; r < s; r++)
                temp[r] = J[r];
        }
    }
    //reduce err
    reduce_err(&err, barrier);

    if (J)
        delete[] J;
    return err;
}





int fill_MSR(double *A, unsigned int *I, unsigned int nx, unsigned int ny, 
             int p, int t, pthread_barrier_t *barrier)
{
    unsigned int size = MSR_offD_size(nx, ny);
    unsigned int n = (nx + 1) * (ny + 1);
    // unsigned int i, j;
    unsigned int l;
    unsigned int s;
    unsigned int *J, *temp;
    double *Am, *atemp;
    unsigned int r = n + 1;
    int err = 0;


    J = new unsigned int[offDMax(nx, ny)];
    if (!J)
        err = -1;
    Am = new double[offDMax(nx, ny)];
    if (!Am)
        err = -2;
    //reduce err
    reduce_err(&err, barrier);

    if (!err)
    {
        if (t == 0)
        {
            for (l = 0; l < n; l++)
            {
                s = offDNum(nx, ny, l);
                I[l] = r;
                r += s;
            }
            I[l] = r;
            // printf("AAAAA %d %d %d %d\n", I[n], I[n - 1], n, size + n);
            if(r != n + size + 1)
                err = -3;
        }
        //reduce err 
        reduce_err(&err, barrier);
    }

    if (!err)
    {
        for (l = tRowsBeg(n, p, t); l < tRowsEnd(n, p, t); l++)
        {
            // printf("%d : %d\n", l, I[n]);
            s = MSR_offD_fill(nx, ny, l, J, A[l], Am);
            if (s != I[l + 1] - I[l])
            {
                // printf("l = %d, s = %d I[l] = %d I[l + 1] = %d\n", l, s, I[l], I[l + 1]);
                err = -4;
            }
            temp = I + I[l];
            atemp = A + I[l];
            for (unsigned int r = 0; r < s; r++)
            {
                temp[r] = J[r];
                atemp[r] = Am[r];
            }
        }
    }
    //reduce err
    reduce_err(&err, barrier);


    if (J)
        delete[] J;
    if (Am)
        delete[] Am;
    return err;
}






int fill_pts(unsigned int nx, unsigned int ny, double x1, double x2, double y1, double y2,
             unsigned int i, unsigned int j, double (*f)(double, double), double fMax, int outlier, double *pts)
{
    double xm = 0., xp = 0., ym = 0., yp = 0., x0 = 0., y0 = 0.;
    double delta_x = (x2 - x1) / nx, delta_y = (y2 - y1) / ny;
    int err = 0;
    int spy = -1;

    if ((i > 1) && (i < nx - 1))
    {
        x0 = x1 + i * delta_x;
        xm = x1 + (i - 1) * delta_x;
        xp = x1 + (i + 1) * delta_x;
    }
    else if (i == 0)
    {
        x0 = x1;
        xm = x1;
        xp = x1 + delta_x;
    }
    else if (i == 1)
    {
        x0 = x1 + delta_x;
        xm = x1;
        xp = x1 + 2 * delta_x;
    }
    else if (i == nx)
    {
        x0 = x2;
        xm = x1 + (nx - 1) * delta_x;
        xp = x2;
    }
    else if (i == nx - 1)
    {
        x0 = x1 + (nx - 1) * delta_x;
        xm = x1 + (nx - 2) * delta_x;
        xp = x2;
    }
    else
    {
        err = -1;
    }


    if ((j > 1) && (j < ny - 1))
    {
        y0 = y1 + j * delta_y;
        ym = y1 + (j - 1) * delta_y;
        yp = y1 + (j + 1) * delta_y;
    }
    else if (j == 0)
    {
        y0 = y1;
        ym = y1;
        yp = y1 + delta_y;
    }
    else if (j == 1)
    {
        y0 = y1 + delta_y;
        ym = y1;
        yp = y1 + 2 * delta_y;
    }
    else if (j == ny)
    {
        y0 = y2;
        ym = y1 + (ny - 1) * delta_y;
        yp = y2;
    }
    else if (j == ny - 1)
    {
        y0 = y1 + (ny - 1) * delta_y;
        ym = y1 + (ny - 2) * delta_y;
        yp = y2;
    }
    else
    {
        err = -2;
    }

    // printf("(%d, %d) : (%e, %e) :: %e\n", i, j, x0, y0, f(x0, y0));

    pts[0] = f(x0, y0);
    pts[1] = f(xp, y0);
    pts[2] = f(x0, ym);
    pts[3] = f(xm, ym);
    pts[4] = f(xm, y0);
    pts[5] = f(x0, yp);
    pts[6] = f(xp, yp);
    pts[7] = f((x0 + xp) / 2., y0);
    pts[8] = f(x0, (y0 + ym) / 2.);
    pts[9] = f((x0 + xp) / 2., (y0 + ym) / 2.);
    pts[10] = f((x0 + xm) / 2., (y0 + ym) / 2.);
    pts[11] = f((x0 + xm) / 2., ym);
    pts[12] = f((x0 + xm) / 2., y0);
    pts[13] = f(xm, (y0 + ym) / 2.);
    pts[14] = f(x0, (y0 + yp) / 2.);
    pts[15] = f((x0 + xm) / 2., (y0 + yp) / 2.);
    pts[16] = f((x0 + xp) / 2., (y0 + yp) / 2.);
    pts[17] = f((x0 + xp) / 2., yp);
    pts[18] = f(xp, (y0 + yp) / 2.);


    if (outlier != 0)
    {
        if ((i == (nx + 1) / 2) && (j == (ny + 1) / 2))
            spy = 0;
        if (((i + 1) == (nx + 1) / 2) && (j == (ny + 1) / 2))
            spy = 1;
        if ((i == (nx + 1) / 2) && ((j - 1) == (ny + 1) / 2))
            spy = 2;
        if (((i - 1) == (nx + 1) / 2) && ((j - 1) == (ny + 1) / 2))
            spy = 3;
        if (((i - 1) == (nx + 1) / 2) && (j == (ny + 1) / 2))
            spy = 4;
        if ((i == (nx + 1) / 2) && ((j + 1) == (ny + 1) / 2))
            spy = 5;
        if (((i + 1) == (nx + 1) / 2) && ((j + 1) == (ny + 1) / 2))
            spy = 6;

        if (spy != -1)
            pts[spy] += fMax * outlier;
    }

    return err;
}


double B_row(unsigned int nx, unsigned int ny, double x1, double x2, double y1, double y2,
             unsigned int l, double (*f)(double, double), double fMax, int outlier, int &err)
{
    unsigned int i, j;
    double ret;
    double pts[19];

    l2ij(nx, ny, i, j, l);
    err = fill_pts(nx, ny, x1, x2, y1, y2, i, j, f, fMax, outlier, pts);
    if (err)
        printf("B_row error %d, %d, %d, %d\n", err, l, i, j);
    else
    {
        if ( ((i > 0) && (i < nx)) && ((j > 0) && (j < ny)) )
        {
            ret = 36. * pts[0] + 20. * ( (pts[7] + pts[8]) + (pts[10] + pts[12]) + (pts[14] + pts[16]) ) + 
                4. * ( (pts[9] + pts[11]) + (pts[13] + pts[15]) + (pts[17] + pts[18]) ) + 
                2. * ( (pts[1] + pts[2]) + (pts[3] + pts[4]) + (pts[5] + pts[6]) );
            return ret / 192.;
        }


        if ( (i == 0) && ((j > 0) && (j < ny)) )
        {
            ret = 18. * pts[0] + 20. * (pts[7] + pts[16]) + 10. * (pts[8] + pts[14]) + 
                4. * (pts[9] + pts[17] + pts[18]) + 2. * (pts[1] + pts[6]) + (pts[2] + pts[5]);
            return ret / 192.;
        }

        if ( (i == nx) && ((j > 0) && (j < ny)) )
        {
            ret = 18. * pts[0] + 20. * (pts[10] + pts[12]) + 10. * (pts[8] + pts[14]) + 
                4. * (pts[11] + pts[13] + pts[15]) + 2. * (pts[3] + pts[4]) + (pts[2] + pts[5]);
            return ret / 192.;
        }

        if ( (j == 0) && ((i > 0) && (i < nx)) )
        {
            ret = 18. * pts[0] + 20. * (pts[14] + pts[16]) + 10. * (pts[7] + pts[12]) + 
                4. * (pts[15] + pts[17] + pts[18]) + 2. * (pts[5] + pts[6]) + (pts[1] + pts[4]);
            return ret / 192.;
        }

        if ( (j == ny) && ((i > 0) && (i < nx)) )
        {
            ret = 18. * pts[0] + 20. * (pts[8] + pts[10]) + 10. * (pts[7] + pts[12]) + 
                4. * (pts[9] + pts[11] + pts[13]) + 2. * (pts[2] + pts[3]) + (pts[1] + pts[4]);
            return ret / 192.;
        }


        if ( (i == 0) && (j == 0) )
        {
            ret = 12. * pts[0] + 20. * pts[16] + 10. * (pts[7] + pts[14]) + 
                4. * (pts[17] + pts[18]) + 2. * pts[6] + (pts[1] + pts[5]);
            return ret / 192.;
        }

        if ( (i == nx) && (j == ny) )
        {
            ret = 12. * pts[0] + 20. * pts[10] + 10. * (pts[8] + pts[12]) + 
                4. * (pts[11] + pts[13]) + 2. * pts[3] + (pts[2] + pts[4]);
            return ret / 192.;
        }


        if ( (i == 0) && (j == ny) )
        {
            ret = 6. * pts[0] + 10. * (pts[7] + pts[8]) + 4. * pts[9] + (pts[1] + pts[2]);
            return ret / 192.;
        }

        if ( (i == nx) && (j == 0) )
        {
            ret = 6. * pts[0] + 10. * (pts[12] + pts[14]) + 4. * pts[15] + (pts[4] + pts[5]);
            return ret / 192.;
        }
    }


    return -1.;
}


int fill_B(double *B, unsigned int nx, unsigned int ny, double x1, double x2, double y1, double y2,
           double (*f)(double, double), double fMax, int outlier, int p, int t, pthread_barrier_t *barrier)
{
    unsigned int n = (nx + 1) * (ny + 1);
    unsigned int l;
    double temp;
    int err = 0, ret = 0;

    for (l = tRowsBeg(n, p, t); l < tRowsEnd(n, p, t); l++)
    {
        temp = B_row(nx, ny, x1, x2, y1, y2, l, f, fMax, outlier, err);
        if (err == 0)
            B[l] = temp;
        else
            ret = -1;
    }
    //reduce err
    reduce_err(&ret, barrier);
    
    return ret;
}








double B_row_debug(unsigned int nx, unsigned int ny, double x1, double x2, double y1, double y2,
    unsigned int l, double (*f)(double, double), double fMax, int outlier, int &err)
{
    unsigned int i, j;
    double ret;
    double pts[19];

    l2ij(nx, ny, i, j, l);
    err = fill_pts(nx, ny, x1, x2, y1, y2, i, j, f, fMax, outlier, pts);
    if (err)
        printf("B_row error\n");
    else
    {
        if ( ((i > 0) && (i < nx)) && ((j > 0) && (j < ny)) )
        {
            ret = 0.5 * pts[0] + (1. / 12.) * ( (pts[1] + pts[2]) + (pts[3] + pts[4]) + (pts[5] + pts[6]) );
            // for (int o = 0; o < 7; o++)
                // printf("center (%d, %d) %d : %e\n", i, j, o, pts[o]);
            return ret;
        }


        if ( (i == 0) && ((j > 0) && (j < ny)) )
        {
            ret = (1. / 4.) * pts[0] + (1. / 12.) * ( pts[1] + pts[6] ) + (1. / 24.) * (pts[2] + pts[5]);
            return ret;
        }

        if ( (i == nx) && ((j > 0) && (j < ny)) )
        {
            ret = (1. / 4.) * pts[0] + (1. / 12.) * ( pts[3] + pts[4] ) + (1. / 24.) * (pts[2] + pts[5]);
            return ret;
        }

        if ( (j == 0) && ((i > 0) && (i < nx)) )
        {
            ret = (1. / 4.) * pts[0] + (1. / 12.) * ( pts[5] + pts[6] ) + (1. / 24.) * (pts[1] + pts[4]);
            // printf("hhhhhhhhhhhhhhhhhhh %e, %e, %e, %e, %e\n", ret, pts[0], pts[4]);
            return ret;
        }

        if ( (j == ny) && ((i > 0) && (i < nx)) )
        {
            ret = (1. / 4.) * pts[0] + (1. / 12.) * ( pts[2] + pts[3] ) + (1. / 24.) * (pts[1] + pts[4]);
            return ret;
        }


        if ( (i == 0) && (j == 0) )
        {
            ret = (1. / 6.) * pts[0] + (1. / 12.) * (pts[6]) + (1. / 24.) * (pts[1] + pts[5]);
            return ret;
        }

        if ( (i == nx) && (j == ny) )
        {
            ret = (1. / 6.) * pts[0] + (1. / 12.) * (pts[3]) + (1. / 24.) * (pts[2] + pts[4]);
            return ret;
        }


        if ( (i == 0) && (j == ny) )
        {
            ret = (1. / 12.) * pts[0] + (1. / 24.) * (pts[1] + pts[2]);
            return ret;
        }

        if ( (i == nx) && (j == 0) )
        {
            ret = (1. / 12.) * pts[0] + (1. / 24.) * (pts[4] + pts[5]);
            return ret;
        }
    }


    return -1.;
}