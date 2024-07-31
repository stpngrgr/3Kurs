#ifndef __APPROX_SURFACE_H__
#define __APPROX_SURFACE_H__

#define W_RES_DEFAULT 10
#define B_RES_DEFAULT 10


class renderSurface
{
    friend class approxWindow;

    private:
        double (*f)(double);
        float *vertexArray;
        unsigned int *indexArray;
        unsigned int wRes = W_RES_DEFAULT;
        unsigned int bRes = B_RES_DEFAULT;
        unsigned int err = 0;

    public:
        renderSurface(double x1, double x2, double y1, double y2, double(*func)(double), unsigned int *reterr);

    private:
        int update(double x1, double x2, double y1, double y2, double(*func)(double), unsigned int *reterr);
        ~renderSurface();

};

#endif