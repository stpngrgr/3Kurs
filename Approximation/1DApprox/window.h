#ifndef __WINDOW_H__
#define __WINDOW_H__

#include <QtWidgets/QtWidgets>
#include <QtWidgets/QShortcut>
#include <QKeySequence>
#include <QString>
#include "approx.h"

class ApproxWindow : public QWidget
{
    Q_OBJECT

    private:
        int funcID, graphID;
        const char *funcName;
        double (*f)(double);
        double a, b;
        int s, p;
        unsigned int n, max_n, drawPt, incr_flg;
        double F_Max = 0.;
        double *X1, *F1, *A1;
        double *X2, *F2, *A2, *Amen;
        QShortcut* hotkeys[8];

        // int Ffill(); 

    public:
        ApproxWindow(QWidget *parent);
        int parse_command_line(int argc, char **argv);
        QSize minimumSizeHint() const;
        QSize sizeHint() const;
        // int initialize(double a, double b, int n);
        QString getInfoString(double F_m, double F_M);
        double firstApproxFunc(double x);
        double secondApproxFunc(double x);
        double func(double x);

    private:
        int firstApproxSetup();
        int secondApproxSetup();
        int approxMemFree();
        int updApprox();
        ~ApproxWindow();

    public slots:
        void func_cycle();
        void graph_cycle();
        void s_inc();
        void s_decr();
        void n_inc();
        void n_decr();
        void p_inc();
        void p_decr();
        void close();

    protected:
        void paintEvent(QPaintEvent *event);

};

#endif