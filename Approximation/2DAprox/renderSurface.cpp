#include "window.h"


approxWindow::renderSurface::renderSurface(approxWindow* _owner)
{
    double *vptr;
    unsigned int *iptr;
    double delta_x, delta_y;
    double x1 = 0., x2 = 0., y1 = 0., y2 = 0.;

    err = 0;

    this->owner = _owner;
    if (!owner)
    {
        err = err | 1;
    }
    else
    {
        x1 = owner->x1;
        x2 = owner->x2;
        y1 = owner->y1;
        y2 = owner->y2;
    }
    vertexArray = new double[wRes * bRes * 3];
    if (!vertexArray)
    {
        err = err | 2;
    }
    indexArray = new unsigned int[(wRes - 1) * (bRes - 1) * 6];
    if (!indexArray)
    {
        err = err | 4;
    }
    if ((x2 - x1 < 1.e-6) || (y2 - y1 < 1.e-6))
    {
        err = err | 8;
    }
    if (!err)
    {   
        {   //vertexArray initialisation
            delta_x = (x2 - x1) / (wRes - 1);
            delta_y = (y2 - y1) / (bRes - 1);
            zMin = zMax = 0.;
            
            vptr = vertexArray;
            vptr[0] = x1;
            vptr[2] = y1;
            vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
            if (vptr[1] < zMin)
                zMin = vptr[1];
            if (vptr[1] > zMax)
                zMax = vptr[1];
            for (unsigned int i = 1; i < wRes - 1; i++)
            {
                vptr += 3;
                vptr[0] = x1 + delta_x * i;
                vptr[2] = y1;
                vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                if (vptr[1] < zMin)
                    zMin = vptr[1];
                if (vptr[1] > zMax)
                    zMax = vptr[1];
            }
            vptr += 3;
            vptr[0] = x2;
            vptr[2] = y1;
            vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
            if (vptr[1] < zMin)
                zMin = vptr[1];
            if (vptr[1] > zMax)
                zMax = vptr[1];

            for (unsigned int j = 1; j < bRes - 1; j++)
            {
                vptr += 3;
                vptr[0] = x1;
                vptr[2] = y1 + delta_y * j;
                vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                if (vptr[1] < zMin)
                    zMin = vptr[1];
                if (vptr[1] > zMax)
                    zMax = vptr[1];
                for (unsigned int i = 1; i < wRes - 1; i++)
                {
                    vptr += 3;
                    vptr[0] = x1 + delta_x * i;
                    vptr[2] = y1 + delta_y * j;
                    vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                    if (vptr[1] < zMin)
                        zMin = vptr[1];
                    if (vptr[1] > zMax)
                        zMax = vptr[1];
                }
                vptr += 3;
                vptr[0] = x2;
                vptr[2] = y1 + delta_y * j;
                vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                if (vptr[1] < zMin)
                    zMin = vptr[1];
                if (vptr[1] > zMax)
                    zMax = vptr[1];
            }

            vptr += 3;
            vptr[0] = x1;
            vptr[2] = y2;
            vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
            if (vptr[1] < zMin)
                zMin = vptr[1];
            if (vptr[1] > zMax)
                zMax = vptr[1];
            for (unsigned int i = 1; i < wRes - 1; i++)
            {
                vptr += 3;
                vptr[0] = x1 + delta_x * i;
                vptr[2] = y2;
                vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                if (vptr[1] < zMin)
                    zMin = vptr[1];
                if (vptr[1] > zMax)
                    zMax = vptr[1];
            }
            vptr += 3;
            vptr[0] = x2;
            vptr[2] = y2;
            vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
            if (vptr[1] < zMin)
                zMin = vptr[1];
            if (vptr[1] > zMax)
                zMax = vptr[1];
        }


        {   //indexArray initialisation
            iptr = indexArray;

            for (unsigned int j = 0; j < bRes - 1; j++)
            {
                iptr[0] = j * wRes;
                iptr[1] = (j + 1) * wRes;
                iptr[2] = iptr[1] + 1;
                iptr[3] = iptr[0];
                iptr[4] = iptr[2];
                iptr[5] = iptr[0] + 1;
                iptr += 6; 
                for (unsigned int i = 1; i < wRes - 1; i++)
                {
                    iptr[0] = iptr[-1];
                    iptr[1] = iptr[-2];
                    iptr[2] = iptr[1] + 1;
                    iptr[3] = iptr[0];
                    iptr[4] = iptr[2];
                    iptr[5] = iptr[0] + 1;
                    iptr += 6; 
                }
            }
        }
    }
}


unsigned int approxWindow::renderSurface::update()
{
    double *vptr;
    double delta_x, delta_y;
    double x1, x2, y1, y2;

    x1 = owner->x1;
    x2 = owner->x2;
    y1 = owner->y1;
    y2 = owner->y2;

    if ((x2 - x1 >= 1.e-6) && (y2 - y1 >= 1.e-6))
    {   
        {   //vertexArray initialisation
            delta_x = (x2 - x1) / (wRes - 1);
            delta_y = (y2 - y1) / (bRes - 1);
            zMin = zMax = 0;

            vptr = vertexArray;
            vptr[0] = x1;
            vptr[2] = y1;
            vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
            if (vptr[1] < zMin)
                zMin = vptr[1];
            if (vptr[1] > zMax)
                zMax = vptr[1];
            for (unsigned int i = 1; i < wRes - 1; i++)
            {
                vptr += 3;
                vptr[0] = x1 + delta_x * i;
                vptr[2] = y1;
                vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                if (vptr[1] < zMin)
                    zMin = vptr[1];
                if (vptr[1] > zMax)
                    zMax = vptr[1];
            }
            vptr += 3;
            vptr[0] = x2;
            vptr[2] = y1;
            vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
            if (vptr[1] < zMin)
                zMin = vptr[1];
            if (vptr[1] > zMax)
                zMax = vptr[1];

            for (unsigned int j = 1; j < bRes - 1; j++)
            {
                vptr += 3;
                vptr[0] = x1;
                vptr[2] = y1 + delta_y * j;
                vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                if (vptr[1] < zMin)
                    zMin = vptr[1];
                if (vptr[1] > zMax)
                    zMax = vptr[1];
                for (unsigned int i = 1; i < wRes - 1; i++)
                {
                    vptr += 3;
                    vptr[0] = x1 + delta_x * i;
                    vptr[2] = y1 + delta_y * j;
                    vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                    if (vptr[1] < zMin)
                        zMin = vptr[1];
                    if (vptr[1] > zMax)
                        zMax = vptr[1];
                }
                vptr += 3;
                vptr[0] = x2;
                vptr[2] = y1 + delta_y * j;
                vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                if (vptr[1] < zMin)
                    zMin = vptr[1];
                if (vptr[1] > zMax)
                    zMax = vptr[1];
            }

            vptr += 3;
            vptr[0] = x1;
            vptr[2] = y2;
            vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
            if (vptr[1] < zMin)
                zMin = vptr[1];
            if (vptr[1] > zMax)
                zMax = vptr[1];
            for (unsigned int i = 1; i < wRes - 1; i++)
            {
                vptr += 3;
                vptr[0] = x1 + delta_x * i;
                vptr[2] = y2;
                vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                if (vptr[1] < zMin)
                    zMin = vptr[1];
                if (vptr[1] > zMax)
                    zMax = vptr[1];
            }
            vptr += 3;
            vptr[0] = x2;
            vptr[2] = y2;
            vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
            if (vptr[1] < zMin)
                zMin = vptr[1];
            if (vptr[1] > zMax)
                zMax = vptr[1];
        }
    }
    else
        return -1;

    return 0;
}   


unsigned int approxWindow::renderSurface::vertexCount()
{
    return wRes * bRes * 3;
}


unsigned int approxWindow::renderSurface::indexCount()
{
    return (wRes - 1) * (bRes - 1) * 6;
}




int approxWindow::renderSurface::update_fBounds(double &fMin, double &fMax, double (*f)(double, double))
{
    double vptr[3];
    double delta_x, delta_y;
    double x1, x2, y1, y2;

    x1 = owner->x1;
    x2 = owner->x2;
    y1 = owner->y1;
    y2 = owner->y2;

    if ((x2 - x1 >= 1.e-6) && (y2 - y1 >= 1.e-6))
    {   
        {   //vertexArray initialisation
            delta_x = (x2 - x1) / (wRes - 1);
            delta_y = (y2 - y1) / (bRes - 1);
            fMin = fMax = 0.;

            vptr[0] = x1;
            vptr[2] = y1;
            vptr[1] = f(vptr[0], vptr[2]);
            if (vptr[1] < fMin)
                fMin = vptr[1];
            if (vptr[1] > fMax)
                fMax = vptr[1];
            for (unsigned int i = 1; i < wRes - 1; i++)
            {
                vptr[0] = x1 + delta_x * i;
                vptr[2] = y1;
                vptr[1] = f(vptr[0], vptr[2]);
                if (vptr[1] < fMin)
                    fMin = vptr[1];
                if (vptr[1] > fMax)
                    fMax = vptr[1];
            }
            vptr[0] = x2;
            vptr[2] = y1;
            vptr[1] = f(vptr[0], vptr[2]);
            if (vptr[1] < fMin)
                fMin = vptr[1];
            if (vptr[1] > fMax)
                fMax = vptr[1];

            for (unsigned int j = 1; j < bRes - 1; j++)
            {
                vptr[0] = x1;
                vptr[2] = y1 + delta_y * j;
                vptr[1] = f(vptr[0], vptr[2]);
                if (vptr[1] < fMin)
                    fMin = vptr[1];
                if (vptr[1] > fMax)
                    fMax = vptr[1];
                for (unsigned int i = 1; i < wRes - 1; i++)
                {
                    vptr[0] = x1 + delta_x * i;
                    vptr[2] = y1 + delta_y * j;
                    vptr[1] = owner->renderFunc(vptr[0], vptr[2]);
                    if (vptr[1] < fMin)
                        fMin = vptr[1];
                    if (vptr[1] > fMax)
                        fMax = vptr[1];
                }
                vptr[0] = x2;
                vptr[2] = y1 + delta_y * j;
                vptr[1] = f(vptr[0], vptr[2]);
                if (vptr[1] < fMin)
                    fMin = vptr[1];
                if (vptr[1] > fMax)
                    fMax = vptr[1];
            }

            vptr[0] = x1;
            vptr[2] = y2;
            vptr[1] = f(vptr[0], vptr[2]);
            if (vptr[1] < fMin)
                fMin = vptr[1];
            if (vptr[1] > fMax)
                fMax = vptr[1];
            for (unsigned int i = 1; i < wRes - 1; i++)
            {
                vptr[0] = x1 + delta_x * i;
                vptr[2] = y2;
                vptr[1] = f(vptr[0], vptr[2]);
                if (vptr[1] < fMin)
                    fMin = vptr[1];
                if (vptr[1] > fMax)
                    fMax = vptr[1];
            }
            vptr[0] = x2;
            vptr[2] = y2;
            vptr[1] = f(vptr[0], vptr[2]);
            if (vptr[1] < fMin)
                fMin = vptr[1];
            if (vptr[1] > fMax)
                fMax = vptr[1];
        }
    }
    else
        return -1;

    return 0;
}
                




int approxWindow::renderSurface::update_fBounds(double &fMin, double &fMax)
{
    double vptr[8];
    double delta_x, delta_y;
    double x1, x2, y1, y2;
    double temp1 = 0., temp2 = 0.;
    double yavrg1 = 0., yavrg2 = 0.;

    x1 = owner->x1;
    x2 = owner->x2;
    y1 = owner->y1;
    y2 = owner->y2;

    if ((x2 - x1 >= 1.e-6) && (y2 - y1 >= 1.e-6))
    {   
        {   //vertexArray initialisation
            delta_x = (x2 - x1) / (wRes - 1);
            delta_y = (y2 - y1) / (bRes - 1);
            fMin = fMax = 0.;

            vptr[0] = x1;
            vptr[1] = y1;
            vptr[2] = x1;
            vptr[3] = y1 + delta_y;
            vptr[4] = x1 + delta_x;
            vptr[5] = y1 + delta_y;
            vptr[6] = x1 + delta_x;
            vptr[7] = y1;
            yavrg1 = (vptr[1] + vptr[3] + vptr[5]) / 3.;
            yavrg2 = (vptr[1] + vptr[5] + vptr[7]) / 3.;
            temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
            temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
            if (temp1 < fMin)
                fMin = temp1;
            if (temp1 > fMax)
                fMax = temp1;
            if (temp2 < fMin)
                fMin = temp2;
            if (temp2 > fMax)
                fMax = temp2;
            for (unsigned int i = 2; i < wRes - 1; i++)
            {
                vptr[0] = vptr[6];
                vptr[2] = vptr[4];
                vptr[4] = x1 + i * delta_x;
                vptr[6] = x1 + i * delta_x;
                temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
                temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
                if (temp1 < fMin)
                    fMin = temp1;
                if (temp1 > fMax)
                    fMax = temp1;
                if (temp2 < fMin)
                    fMin = temp2;
                if (temp2 > fMax)
                    fMax = temp2;
            }
            vptr[0] = x1 + (wRes - 2) * delta_x;
            vptr[2] = x1 + (wRes - 2) * delta_x;
            vptr[4] = x2;
            vptr[6] = x2;
            temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
            temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
            if (temp1 < fMin)
                fMin = temp1;
            if (temp1 > fMax)
                fMax = temp1;
            if (temp2 < fMin)
                fMin = temp2;
            if (temp2 > fMax)
                fMax = temp2;

            for (unsigned int j = 2; j < bRes - 1; j++)
            {
                vptr[0] = x1;
                vptr[1] = y1 + (j - 1) * delta_y;
                vptr[2] = x1;
                vptr[3] = y1 + j * delta_y;
                vptr[4] = x1 + delta_x;
                vptr[5] = y1 + j * delta_y;
                vptr[6] = x1 + delta_x;
                vptr[7] = y1 + (j - 1) * delta_y;
                yavrg1 = (vptr[1] + vptr[3] + vptr[5]) / 3.;
                yavrg2 = (vptr[1] + vptr[5] + vptr[7]) / 3.;
                temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
                temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
                if (temp1 < fMin)
                    fMin = temp1;
                if (temp1 > fMax)
                    fMax = temp1;
                if (temp2 < fMin)
                    fMin = temp2;
                if (temp2 > fMax)
                    fMax = temp2;
                for (unsigned int i = 2; i < wRes - 1; i++)
                {
                    vptr[0] = vptr[6];
                    vptr[2] = vptr[4];
                    vptr[4] = x1 + i * delta_x;
                    vptr[6] = x1 + i * delta_x;
                    temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
                    temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
                    if (temp1 < fMin)
                        fMin = temp1;
                    if (temp1 > fMax)
                        fMax = temp1;
                    if (temp2 < fMin)
                        fMin = temp2;
                    if (temp2 > fMax)
                        fMax = temp2;
                }
                vptr[0] = x1 + (wRes - 2) * delta_x;
                vptr[2] = x1 + (wRes - 2) * delta_x;
                vptr[4] = x2;
                vptr[6] = x2;
                temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
                temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
                if (temp1 < fMin)
                    fMin = temp1;
                if (temp1 > fMax)
                    fMax = temp1;
                if (temp2 < fMin)
                    fMin = temp2;
                if (temp2 > fMax)
                    fMax = temp2;
            }

            vptr[0] = x1;
            vptr[1] = y1 + (bRes - 2) * delta_y;
            vptr[2] = x1;
            vptr[3] = y2;
            vptr[4] = x1 + delta_x;
            vptr[5] = y2;
            vptr[6] = x1 + delta_x;
            vptr[7] = y1 + (bRes - 2) * delta_y;
            yavrg1 = (vptr[1] + vptr[3] + vptr[5]) / 3.;
            yavrg2 = (vptr[1] + vptr[5] + vptr[7]) / 3.;
            temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
            temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
            if (temp1 < fMin)
                fMin = temp1;
            if (temp1 > fMax)
                fMax = temp1;
            if (temp2 < fMin)
                fMin = temp2;
            if (temp2 > fMax)
                fMax = temp2;
            for (unsigned int i = 2; i < wRes - 1; i++)
            {
                vptr[0] = vptr[6];
                vptr[2] = vptr[4];
                vptr[4] = x1 + i * delta_x;
                vptr[6] = x1 + i * delta_x;
                temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
                temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
                if (temp1 < fMin)
                    fMin = temp1;
                if (temp1 > fMax)
                    fMax = temp1;
                if (temp2 < fMin)
                    fMin = temp2;
                if (temp2 > fMax)
                    fMax = temp2;
            }
            vptr[0] = x1 + (wRes - 2) * delta_x;
            vptr[2] = x1 + (wRes - 2) * delta_x;
            vptr[4] = x2;
            vptr[6] = x2;
            temp1 = owner->renderFunc((vptr[0] + vptr[2] + vptr[4]) / 3., yavrg1);
            temp2 = owner->renderFunc((vptr[0] + vptr[4] + vptr[6]) / 3., yavrg2);
            if (temp1 < fMin)
                fMin = temp1;
            if (temp1 > fMax)
                fMax = temp1;
            if (temp2 < fMin)
                fMin = temp2;
            if (temp2 > fMax)
                fMax = temp2;

        }
    }
    else
        return -1;

    return 0;
}







approxWindow::renderSurface::~renderSurface()
{
    if (!(err & 2))
    {
        delete[] vertexArray;
    }
    if (!(err & 4))
    {
        delete[] indexArray;
    }
}










// approxWindow::renderSurface::renderSurface(approxWindow* _owner)
// {


//     this->owner = _owner;
//     if (!owner)
//     {
//         err | 1;
//     }
//     vertexArray = new double[9];
//     if (!vertexArray)
//     {
//         err | 2;
//     }
//     indexArray = new unsigned int[3];
//     if (!indexArray)
//     {
//         err | 4;
//     }
//     if (!err)
//     {
//         indexArray[0] = 0;
//         indexArray[1] = 1;
//         indexArray[2] = 2;
//         vertexArray[0] = -0.5f;
//         vertexArray[1] = -0.5f;
//         vertexArray[2] = -2.0f;
//         vertexArray[3] = 0.5f;
//         vertexArray[4] = -0.5f;
//         vertexArray[5] = -2.0f;
//         vertexArray[6] = 0.0f;
//         vertexArray[7] = 0.5f;
//         vertexArray[8] = -30.0f;
//     }
// }


// unsigned int approxWindow::renderSurface::vertexCount()
// {
//     return 9;
// }


// unsigned int approxWindow::renderSurface::indexCount()
// {
//     return 3;
// }