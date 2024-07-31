#ifndef __THREAD_MANAGER_H__
#define __THREAD_MANAGER_H__
#include <QtWidgets/QApplication>
#include <QObject>
#include <QEvent>


#include <unistd.h>
#include <pthread.h>

// #include "window.h"
#include "Solve.h"
#include "fill_matrix.h"




struct thread_args
{
    class calcThreadManager *manager;
    class guiThreadComm *comm;
    int t;
    int p;
};




void* thread_do(void* args);


void* thread_main_loop(void* args);






class timeUpdEvent : public QEvent
{

    public:
        timeUpdEvent(QEvent::Type type);
        ~timeUpdEvent();
};





class threadManager
{
    public:
        const QEvent::Type timeUpdEventType;

    public:
        threadManager();
        int wakeGUI();
        int attachApp(QObject *app);
        ~threadManager();

    private:
       QObject *app = nullptr;
};





void* thread_main_loop(void* args);


class calcThreadManager
{
    // friend class approxWindow;

    public:
        calcThreadManager(double _x1, double _x2, double _y1, double _y2, int _maxit, double _eps, int _p);
        int getThreadArgs(thread_args* args, guiThreadComm *comm);
        ~calcThreadManager();


    public:
        pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
        pthread_cond_t ready = PTHREAD_COND_INITIALIZER;
        pthread_cond_t rquit = PTHREAD_COND_INITIALIZER;
        pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
        pthread_barrier_t barrier;
        int state = 0; // 0 = no work, 1 = working, -1 = quit
        int count = 0;
        int ready_flg = 0, quit_flg = 0;

    public:
        int p;

        unsigned int nx, ny;
        double x1, x2, y1, y2;

        double *A = nullptr;
        unsigned int *I = nullptr;
        double *b = nullptr, *x = nullptr;
        double *r = nullptr, *u = nullptr, *v = nullptr;
        double* buffer = nullptr;

        double (*f)(double, double);
        double fMax;
        int outlier;

        int err = 0;
        int maxit;
        double eps;
};

class guiThreadComm
{
    public :
        guiThreadComm();
        int wakeGUI();
        int attachApp(QObject *_app);
        ~guiThreadComm();

    public:
        const QEvent::Type wakeEventType;
        QObject *app = nullptr;
        pthread_cond_t wakeupcall = PTHREAD_COND_INITIALIZER;
        pthread_cond_t ready = PTHREAD_COND_INITIALIZER;
        pthread_cond_t esheodin = PTHREAD_COND_INITIALIZER;
        pthread_cond_t rquit = PTHREAD_COND_INITIALIZER;
        pthread_mutex_t guitex = PTHREAD_MUTEX_INITIALIZER;
        int request = 0; // 0 = no request, 1 = recalculation, -1 = quit
        int state = 0; // 0 = ready, 1 = busy, 2 = done, -1 = quitting 
        int err = 0;
        int count = 0;
        int ready_flg = 0, rquit_flg = 0, esheodin_flg = 0;

        double *x = nullptr;
        unsigned int nx, ny;
        double (*f)(double, double);
        double fMax;
        int outlier;
        int iteration;
};

#endif