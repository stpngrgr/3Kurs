#define _USE_MATH_DEFINES

#include <QPainter>
#include <QString>

#include "window.h"
#include "approx.h"
#include <math.h>
#include <stdio.h>

static double f0(double) // f(x) = 1
{
    return 1.;
}


static double f1(double x) // f(x) = x
{
    return x;
}


static double f2(double x) // f(x) = x^2
{
    return x * x;
}


static double f3(double x) // f(x) = x^3
{
    return x * x * x;
}


static double f4(double x) // f(x) = x^4
{
    return x * x * x * x;
}


static double f5(double x) // f(x) = e^x
{
    return exp(x);
}


static double f6(double x) // f(x) = 1/(25x^2 + 1)
{
    return 1 / (25 * x * x + 1);
}


ApproxWindow::ApproxWindow(QWidget *parent) : QWidget (parent)
{
    a = -10.;
    b = 10.;
    drawPt = 1000;
    s = 0;
    p = 0;
    n = 8;
    max_n = n;
    incr_flg = 3;
    f = f1;
    funcName = "k = 1 f(x) = x";
    funcID = 1;
    graphID = 0;

    X1 = F1 = A1 = nullptr;
    X2 = F2 = A2 = Amen = nullptr;

    hotkeys[0] = new QShortcut(QKeySequence("0"), this, SLOT(func_cycle()));
    hotkeys[1] = new QShortcut(QKeySequence("1"), this, SLOT(graph_cycle()));
    hotkeys[2] = new QShortcut(QKeySequence("2"), this, SLOT(s_inc()));
    hotkeys[3] = new QShortcut(QKeySequence("3"), this, SLOT(s_decr()));
    hotkeys[4] = new QShortcut(QKeySequence("4"), this, SLOT(n_inc()));
    hotkeys[5] = new QShortcut(QKeySequence("5"), this, SLOT(n_decr()));
    hotkeys[6] = new QShortcut(QKeySequence("6"), this, SLOT(p_inc()));
    hotkeys[7] = new QShortcut(QKeySequence("7"), this, SLOT(p_decr()));

    updApprox();
    update();
}


int ApproxWindow::parse_command_line(int argc, char **argv)
{
    int err = 0;

    if (argc == 1)
        err = 1;
    else if (argc == 5)
    {
        if ((sscanf(argv[1], "%lf", &a) != 1) || (sscanf(argv[2], "%lf", &b) != 1) || ((b - a) < 1.e-6))
            err = -1;
        if ((sscanf(argv[3], "%d", &n) != 1) || (n < 2))
            err = -2;
        if ((sscanf(argv[4], "%d", &funcID) != 1) || (funcID < 0) || (funcID > 6))
            err = -3;
        else
        {
            funcID =  (funcID != 0 ? funcID - 1 : 6);
            if (n > max_n)
            {
                incr_flg = 3;
                max_n = n;
            }
            func_cycle();
        }
    }
    else
        err = -5;

    return err;
}





void ApproxWindow::func_cycle()
{
    double x, y, delta_x = (b - a) / drawPt;

    funcID = (funcID + 1) % 7;

    switch (funcID)
    {
        case 0:
            f = f0;
            funcName = "k = 0 | f(x) = 1";
            break;
        case 1:
            f = f1;
            funcName = "k = 1 | f(x) = x";
            break;
        case 2:
            f = f2;
            funcName = "k = 2 | f(x) = x^2";
            break;
        case 3:
            f = f3;
            funcName = "k = 3 | f(x) = x^3";
            break;
        case 4:
            f = f4;
            funcName = "k = 4 | f(x) = x^4";
            break;
        case 5:
            f = f5;
            funcName = "k = 5 | f(x) = e^x";
            break;
        case 6:
            f = f6;
            funcName = "k = 6 | f(x) = 1/(25x^2 + 1)";
            break;
        default:
            f = f0;
            funcName = "k = 0 | f(x) = 1";
            perror("func_cycle default case executed\n");
            break;
    }

    F_Max = 0.;
    for (x = a; x - b < 1.e-6; x += delta_x)
    {
       y = fabs(func(x));
        if (y > F_Max)
            F_Max = y;
    }

    updApprox();
    update();
}


void ApproxWindow::graph_cycle()
{
    graphID = (graphID + 1) % 5;
    update();
}


void ApproxWindow::s_inc()
{
    s++;
    update();
}


void ApproxWindow::s_decr()
{
    s--;
    if (s < 0)
        s = 0;
    update();
}


void ApproxWindow::n_inc()
{
    n *= 2;
    if (n > max_n)
    {
        incr_flg = 3;
        max_n = n;
    }
    updApprox();
    update();
}


void ApproxWindow::n_decr()
{
    n /= 2;
    if (n < 2)
        n = 2;
    updApprox();
    update();
}


void ApproxWindow::p_inc()
{
    p++;
    updApprox();
    update();
}


void ApproxWindow::p_decr()
{
    p--;
    updApprox();
    update();
}





QString ApproxWindow::getInfoString(double F_m, double F_M)
{
    QString infoStr = QString("%1 | GraphID = %2 | s = %3 | p = %4 | [F_m, F_M] = [%5, %6] | n = %7").arg(
        funcName, QString::number(graphID), QString::number(s), QString::number(p), QString::number(F_m), QString::number(F_M), 
        QString::number(n));

    return infoStr;
}


double ApproxWindow::firstApproxFunc(double x)
{
    return firstApprox_func(x, a, b, n, X1, A1);
}

double ApproxWindow::secondApproxFunc(double x)
{
    return secondApprox_func(x, a, b, n, X2, A2);
}


double ApproxWindow::func(double x)
{
    return f(x);
}


int ApproxWindow::updApprox()
{
    int err = 0;

    if (n <= 50)
    {
        err = firstApproxSetup();
        if (!err)
            err = firstApprox(n, X1, F1, A1, nullptr);
        if (err)
        {
            printf("First approximation error\n");
            return err;
        }
    }
    err = secondApproxSetup();
    if (!err)
        err = secondApprox(n, X2, F2, A2, Amen);
    if (err)
    {
        printf("Second approximation error\n");
        return err;
    }

    return 0;
}






int ApproxWindow::firstApproxSetup()
{
    if (incr_flg & 1)
    {
        if (X1)
            delete[] (X1);
        X1 = new double[n];
        if (!X1)
        {
            return -1;
        }
        if (F1)
            delete[] (F1);
        F1 = new double[n];
        if (!F1)
        {
            return -1;
        }
        if (A1)
            delete[] (A1);
        A1 = new double[3 * n + 3];
        if (!A1)
        {
            return -1;
        }
        incr_flg = incr_flg & 2;
    }


    for (unsigned int i = 1; i < n + 1; i++)
    {
        X1[i - 1] = ((a + b) + ((b - a) * cos((M_PI * (2 * i - 1)) / (2 * n)))) / 2;
        F1[i - 1] = func(X1[i - 1]);
        if ((p != 0) && ((i - 1) == n / 2))
            F1[i - 1] += p * 0.1 * F_Max;
        // printf ("%lf\n", (M_PI * (2 * i - 1)) / (2 * n));
        // printf("%d : (%10.10lf, %10.10lf)\n", i, X1[i - 1], F1[i - 1]);
    }


    // X1[0] = a;
    // F1[0] = func(a);
    // if ((p != 0) && (0 == n / 2))
    //     F1[0] += p * 0.1 * F_Max;
    // X1[n - 1] = b;
    // F1[n - 1] = func(b);

    // for (int i = 1; i < n - 1; i++)
    // {
    //     // X2[i] = X2[i - 1] + delta_x;
        
    //     X1[i] = ((a + b) + ((b - a) * cos((M_PI * (2 * i - 1)) / (2 * (n - 2))))) / 2;

    //     F1[i] = func(X1[i]);
    //     if ((p != 0) && (i == n / 2))
    //         F1[i] += p * 0.1 * F_Max;
    //     // printf("%d : (%10.10lf, %10.10lf)\n", i, X2[i], F2[i]);
    // }
    
    return 0;
}


int ApproxWindow::secondApproxSetup()
{
    double delta_x = 0.;


    if (incr_flg & 2)
    {
        if (X2)
            delete[] (X2);
        X2 = new double[n];
        if (!X2)
        {
            return -1;
        }
        if (F2)
            delete[] (F2);
        F2 = new double[n];
        if (!F2)
        {
            return -1;
        }
        if (A2)
            delete[] (A2);
        A2 = new double[3 * n + 3];
        if (!A2)
        {
            return -1;
        }
        if (Amen)
            delete[] (Amen);
        Amen = new double[2 * n + 2];
        if (!Amen)
        {
            return -1;
        }
        incr_flg = incr_flg & 1;
    }

    if (n < 2)
        return -1;
    else if (n > 2)
        delta_x = (b - a) / (n - 1);

    X2[0] = a;
    F2[0] = func(a);
    if ((p != 0) && (0 == n / 2))
        F2[0] += p * 0.1 * F_Max;
    X2[n - 1] = b;
    F2[n - 1] = func(b);

    for (unsigned int i = 1; i < n - 1; i++)
    {
        X2[i] = a + delta_x * i;

        // X2[i] = ((a + b) + ((b - a) * cos((M_PI * (2 * i - 1)) / (2 * (n - 2))))) / 2;

        F2[i] = func(X2[i]);
        if ((p != 0) && (i == n / 2))
            F2[i] += p * 0.1 * F_Max;
        // printf("%d : (%10.10lf, %10.10lf)\n", i, X2[i], F2[i]);
    }

    return 0;
}





int ApproxWindow::approxMemFree()
{
    if (X1)
    {
        delete[] (X1);
        X1 = nullptr;
    }
    if (F1)
    {
        delete[] (F1);
        F1 = nullptr;
    }
    if (A1)
    {
        delete[] (A1);
        A1 = nullptr;
    }
    if (X2)
    {
        delete[] (X2);
        X2 = nullptr;
    }
    if (F2)
    {
        delete[] (F2);
        F2 = nullptr;
    }
    if (A2)
    {
        delete[] (A2);
        A2 = nullptr;
    }
    if (Amen)
    {
        delete[] (Amen);
        Amen = nullptr;
    }

    return 0;
}

void ApproxWindow::close()
{
    approxMemFree();
    for (int i = 0; i < 8; i++)
    {
        if (hotkeys[i])
        {
            delete hotkeys[i];
            hotkeys[i] = nullptr;
        }
    }
}





static inline int sc(double x, double a, double b, int h)
{
    double res;
    int ret;
     
    if (b - a > 1.e-16)
    	res = h * (x - a) / (b - a);
    else
        return -1;
        
    ret = (res - ((int)res) <= 0.5 ? (int)res : (int)res + 1);
    if (ret < 0)
    	ret = 0;
    else if (ret > h)
    	ret = h;
    	
    return ret;
}





void ApproxWindow::paintEvent(QPaintEvent * /*event*/)
{
    QPainter painter(this);
    QPen pen_black(Qt::black, 0, Qt::SolidLine);
    QPen pen_blue(Qt::blue, 0, Qt::SolidLine);  
    QPen pen_red(Qt::red, 0, Qt::SolidLine);  
    QPen pen_green(Qt::green, 0, Qt::SolidLine);  
    double as = a / pow(2, (double) s), bs = b / pow(2, (double) s);
    double x1 = as, x2, y1, y2, l1, l2, y11, y22;
    double max_y = 0., min_y = 0.; 
    double F_m = 0., F_M = 0.;
    double delta_x, delta_y;
    int winh = height();
    int winw = width();
    int ptx, pty;

    drawPt = winw;
    delta_x = (bs - as) / winw;

    switch (graphID)
    {
        case 0 :

            for (x2 = as; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y2 = func(x2);
                if (y2 < min_y)
                    min_y = y2;
                if (y2 > max_y)
                    max_y = y2;
            }
            if (n <= 50)
            {
                for (x2 = as; x2 - bs < 1.e-6; x2 += delta_x)
                {
                    y1 = firstApproxFunc(x2);
                    if (y1 < min_y)
                        min_y = y1;
                    if (y1 > max_y)
                        max_y = y1;
                } 
            }

            delta_y = 0.01 * (max_y - min_y);
            F_m = min_y;
            min_y -= delta_y;
            F_M = max_y;
            max_y += delta_y;


            painter.setPen(pen_black);
            pty = winh - sc(0., min_y, max_y, winh);
            painter.drawLine(0, pty, winw, pty);
            ptx = sc(0., as, bs, winw);
            painter.drawLine(ptx, 0, ptx, winh);

            // painter.drawLine(QPointF(-6., 6.), QPointF(-2., 2.));


            painter.setPen(pen_blue);

            y1 = func(x1);
            for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y2 = func(x2);
                painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));
                x1 = x2;
                y1 = y2;
            }
            // x2 = bs;
            // y2 = func(x2);
            // painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));

            if (n <= 50)
            {
                painter.setPen(pen_red);
                x1 = as;
                l1 = firstApproxFunc(x1);
                for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
                {
                    l2 = firstApproxFunc(x2);
                    painter.drawLine(sc(x1, as, bs, winw), winh - sc(l1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(l2, min_y, max_y, winh));
                    x1 = x2;
                    l1 = l2;
                }
                // x2 = bs;
                // l2 = firstApproxFunc(x2);
                // painter.drawLine(sc(x1, as, bs, winw), winh - sc(l1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(l2, min_y, max_y, winh));
            }

            break;
            //case 0 end

        case 1 : 

            for (x2 = as; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y1 = secondApproxFunc(x2);
                y2 = func(x2);
                if (y2 < min_y)
                    min_y = y2;
                if (y2 > max_y)
                    max_y = y2;
                if (y1 < min_y)
                    min_y = y1;
                if (y1 > max_y)
                    max_y = y1;
            }

            delta_y = 0.01 * (max_y - min_y);
            F_m = min_y;
            min_y -= delta_y;
            F_M = max_y;
            max_y += delta_y;


            painter.setPen(pen_black);
            pty = winh - sc(0., min_y, max_y, winh);
            painter.drawLine(0, pty, winw, pty);
            ptx = sc(0., as, bs, winw);
            painter.drawLine(ptx, 0, ptx, winh);

            // painter.drawLine(QPointF(-6., 6.), QPointF(-2., 2.));


            painter.setPen(pen_blue);

            y1 = func(x1);
            for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y2 = func(x2);
                painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));
                x1 = x2;
                y1 = y2;
            }
            // x2 = bs;
            // y2 = func(x2);
            // painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));

            {
                painter.setPen(pen_green);
                x1 = as;
                l1 = secondApproxFunc(x1);
                for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
                {
                    l2 = secondApproxFunc(x2);
                    painter.drawLine(sc(x1, as, bs, winw), winh - sc(l1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(l2, min_y, max_y, winh));
                    x1 = x2;
                    l1 = l2;
                    // printf("\\\\\\\\\\\\\\\\%10.10e, %10.10e\n", x2, l2);
                }
                // x2 = bs;
                // l2 = secondApproxFunc(x2);
                // painter.drawLine(sc(x1, as, bs, winw), winh - sc(l1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(l2, min_y, max_y, winh));
            }

            break;
            //case 1 end

        case 2:


            for (x2 = as; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y11 = secondApproxFunc(x2);
                y2 = func(x2);
                if (y2 < min_y)
                    min_y = y2;
                if (y2 > max_y)
                    max_y = y2;
                if (y11 < min_y)
                    min_y = y11;
                if (y11 > max_y)
                    max_y = y11;
            }
            if (n <= 50)
            {
                for (x2 = as; x2 - bs < 1.e-6; x2 += delta_x)
                {
                    y1 = firstApproxFunc(x2);
                    if (y1 < min_y)
                        min_y = y1;
                    if (y1 > max_y)
                        max_y = y1;
                } 
            }

            delta_y = 0.01 * (max_y - min_y);
            F_m = min_y;
            min_y -= delta_y;
            F_M = max_y;
            max_y += delta_y;


            painter.setPen(pen_black);
            pty = winh - sc(0., min_y, max_y, winh);
            painter.drawLine(0, pty, winw, pty);
            ptx = sc(0., as, bs, winw);
            painter.drawLine(ptx, 0, ptx, winh);

            // painter.drawLine(QPointF(-6., 6.), QPointF(-2., 2.));


            painter.setPen(pen_blue);

            y1 = func(x1);
            for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y2 = func(x2);
                painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));
                x1 = x2;
                y1 = y2;
            }
            // x2 = bs;
            // y2 = func(x2);
            // painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));

            if (n <= 50)
            {
                painter.setPen(pen_red);

                x1 = as;
                y1 = firstApproxFunc(x1);
                for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
                {
                    y2 = firstApproxFunc(x2);
                    painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));
                    x1 = x2;
                    y1 = y2;
                }
                // x2 = bs;
                // y2 = firstApproxFunc(x2);
                // painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));
            }

            painter.setPen(pen_green);

            x1 = as;
            y11 = secondApproxFunc(x1);
            for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y22 = secondApproxFunc(x2);
                painter.drawLine(sc(x1, as, bs, winw), winh - sc(y11, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y22, min_y, max_y, winh));
                x1 = x2;
                y11 = y22;
            }
            // x2 = bs;
            // y22 = secondApproxFunc(x2);
            // painter.drawLine(sc(x1, as, bs, winw), winh - sc(y11, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y22, min_y, max_y, winh));

            break;
            //case 2 end

        case 3:

            // for (x2 = as; x2 - bs < 1.e-6; x2 += delta_x)
            // {
            //     y22 = fabs(secondApproxFunc(x2) - func(x2));
            //     printf("%lf, %lf, %lf\n", secondApproxFunc(x2), func(x2), y22);
            //     if (y22 < min_y)
            //         min_y = y22;
            //     if (y22 > max_y)
            //         max_y = y22;
            // }
            if (n <= 50)
            {
                for (x2 = as; x2 - bs < 1.e-6; x2 += delta_x)
                {
                    y2 = fabs(firstApproxFunc(x2) - func(x2));
                    if (y2 < min_y)
                        min_y = y2;
                    if (y2 > max_y)
                        max_y = y2;
                } 
            }

            delta_y = 0.01 * (max_y - min_y);
            F_m = min_y;
            min_y -= delta_y;
            F_M = max_y;
            max_y += delta_y;


            painter.setPen(pen_black);
            pty = winh - sc(0., min_y, max_y, winh);
            painter.drawLine(0, pty, winw, pty);
            ptx = sc(0., as, bs, winw);
            painter.drawLine(ptx, 0, ptx, winh);

            // painter.drawLine(QPointF(-6., 6.), QPointF(-2., 2.));

            if (n <= 50)
            {
                painter.setPen(pen_red);

                y1 = fabs(firstApproxFunc(x1) - func(x1));
                for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
                {
                    y2 = fabs(firstApproxFunc(x2) - func(x2));
                    painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));
                    x1 = x2;
                    y1 = y2;
                }
                // x2 = bs;
                // y2 = fabs(firstApproxFunc(x2) - func(x2));
                // painter.drawLine(sc(x1, as, bs, winw), winh - sc(y1, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y2, min_y, max_y, winh));
            }


            // painter.setPen(pen_green);

            // x1 = as;
            // y11 = fabs(secondApproxFunc(x1) - func(x1));
            // for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
            // {
            //     y22 = fabs(secondApproxFunc(x2) - func(x2));
            //     painter.drawLine(QPointF(x1, y11), QPointF(x2, y22));
            //     x1 = x2;
            //     y11 = y22;
            // }
            // x2 = bs;
            // y22 = fabs(secondApproxFunc(x2) - func(x2));
            // painter.drawLine(QPointF(x1, y11), QPointF(x2, y22));

            break;
            //case 3 end

        case 4:

            for (x2 = as; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y22 = fabs(secondApproxFunc(x2) - func(x2));
                // printf("%lf, %lf, %lf\n", secondApproxFunc(x2), func(x2), y22);
                if (y22 < min_y)
                    min_y = y22;
                if (y22 > max_y)
                    max_y = y22;
            }

            delta_y = 0.01 * (max_y - min_y);
            F_m = min_y;
            min_y -= delta_y;
            F_M = max_y;
            max_y += delta_y;


            painter.setPen(pen_black);
            pty = winh - sc(0., min_y, max_y, winh);
            painter.drawLine(0, pty, winw, pty);
            ptx = sc(0., as, bs, winw);
            painter.drawLine(ptx, 0, ptx, winh);

            // painter.drawLine(QPointF(-6., 6.), QPointF(-2., 2.));


            painter.setPen(pen_green);

            x1 = as;
            y11 = fabs(secondApproxFunc(x1) - func(x1));
            for (x2 = x1 + delta_x; x2 - bs < 1.e-6; x2 += delta_x)
            {
                y22 = fabs(secondApproxFunc(x2) - func(x2));
                painter.drawLine(sc(x1, as, bs, winw), winh - sc(y11, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y22, min_y, max_y, winh));
                x1 = x2;
                y11 = y22;
            }
            // x2 = bs;
            // y22 = fabs(secondApproxFunc(x2) - func(x2));
            // painter.drawLine(sc(x1, as, bs, winw), winh - sc(y11, min_y, max_y, winh), sc(x2, as, bs, winw), winh - sc(y22, min_y, max_y, winh));

            break;
            //case 3 end

        default : 
            break;
    }


    painter.setPen(pen_blue);
    painter.drawText (0, 20, getInfoString(F_m, F_M));
    printf("Scale : %10.6e\n", (fabs(F_M) > fabs(F_m) ? fabs(F_M) : fabs(F_m)));
}





QSize ApproxWindow::minimumSizeHint () const
{
  return QSize (100, 100);
}


QSize ApproxWindow::sizeHint () const
{
  return QSize (1000, 1000);
}





ApproxWindow::~ApproxWindow()
{
    approxMemFree();
    for (int i = 0; i < 8; i++)
    {
        if (hotkeys[i])
        {
            delete hotkeys[i];
            hotkeys[i] = nullptr;
        }
    }
}
