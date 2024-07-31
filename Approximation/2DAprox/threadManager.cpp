

#include <cstring>
#include <math.h>

#include "threadManager.h"




void* thread_do(void *args)
{
    thread_main_loop(args);
    // printf("exit : %d\n", ((thread_args*)(args))->t);
    return nullptr;
}




// static void display(double *x, unsigned int n)
// {
//     printf("\n\n ---------------- \n\n");
//     for (unsigned int i = 0; i < n; i++)
//     {
//         printf("%e\n", x[i]);
//     }
//     printf("\n\n ---------------- \n\n");
// }



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






int solve(calcThreadManager *manager, unsigned int p, unsigned int t, pthread_barrier_t *barrier)
{
    int err = 0;
    int it;

    double norm = 0., temp = 0.;

    unsigned int nx = manager->nx, ny = manager->ny;
    double x1 = manager->x1, x2 = manager->x2, y1 = manager->y1, y2 = manager->y2;

    double *A = manager->A;
    unsigned int *I = manager->I;
    double *b = manager->b, *x = manager->x;
    double *r = manager->r, *u = manager->u, *v = manager->v;
    double *buffer = manager->buffer;

    double (*f)(double, double) = manager->f;
    double fMax = manager->fMax;
    int outlier = manager->outlier;
    int maxit = manager->maxit;
    double eps = manager->eps;

    unsigned int n = (nx + 1) * (ny + 1);
    double *amt = nullptr;

    amt = new double[(tRowsEnd(n, p, t) - tRowsBeg(n, p, t)) / 2 + 1];
    if (!amt)
    {
        err = -1;
        amt = nullptr;
    }
    //reduce err;
    reduce_err(&err, barrier);

    if (err)
    {
        if (amt)
            delete[] amt;
        amt = nullptr;
        return err;
    }

    err = fill_MSR(A, I, nx, ny, p, t, barrier);
    if (err)
    {
        if (amt)
            delete[] amt;
        amt = nullptr;
        err = -2;
        return err;
    }
    err = fill_B(b, nx, ny, x1, x2, y1, y2, f, fMax, outlier, p, t, barrier);
    if (err)
    {
        if (amt)
            delete[] amt;
        amt = nullptr;
        err = -3;
        return err;
    }



    for (unsigned int i = tRowsBeg(n, p, t); i < tRowsEnd(n, p, t); i++)
    {
        x[i] = 1.;
    }

    if (t == 0)
    {
        for (unsigned int i = 0; i < n; i++)
        {
            temp = A[i];

            for (unsigned int j = I[i]; j < I[i + 1]; j++)
            {
                temp += A[j];
            }

            norm = (fabs(temp) < norm ? norm : fabs(temp));
        }

        for (unsigned int i = 0; i < n; i++)
            if (A[i] < 1.e-16 * norm)
                err = -5;
    }


    reduce_err(&err, barrier);
    if (err)
    {
        if (amt)
            delete[] amt;
        amt = nullptr;
        // printf("error = %d\n", err);
        err = -4;
        return err;
    }



    // display(A, n);
    // display(b, n);

    // printf("%d %e %d \n", n, eps, maxit);

    it = tSolveRestart(A, I, b, x, r, u, v, amt, buffer, n, eps, maxit, p, t, barrier);

    // display(x, n);

    // printf("out %d\n", t);

    delete[] amt;
    amt = nullptr;
    return it;
}





void* thread_main_loop(void* args)
{
    calcThreadManager *manager = ((thread_args*)(args))->manager;
    int t = ((thread_args*)(args))->t, p = ((thread_args*)(args))->p;
    pthread_cond_t *cond = &(manager->cond);
    pthread_cond_t *ready = &(manager->ready);
    pthread_mutex_t *mutex = &(manager->mutex);
    pthread_barrier_t *barrier = &(manager->barrier);

    guiThreadComm *comm = ((thread_args*)(args))->comm;

    int err = 0, ret;
    unsigned int n;
    int quit = 0;

    if (t == 0)
    {


        pthread_mutex_lock(&(comm->guitex));


        while (quit != -1)
        {
            pthread_mutex_lock(mutex);
            manager->count++;
            if (manager->count >= p)
            {
                manager->count = 0;
            }
            else
            {
                // printf("rrthererr\n");
                while (!manager->ready_flg)
                    pthread_cond_wait(ready, mutex);
            }
            manager->ready_flg = 0;
            // printf("every\n");
            pthread_mutex_unlock(mutex);

            comm->ready_flg = 1;
            pthread_cond_broadcast(&(comm->ready));
            pthread_cond_wait(&(comm->wakeupcall), &(comm->guitex));
            comm->ready_flg = 0;
            switch(comm->request)
            {
                case 0: //no request
                    break;

                case -1: //quit
                    quit = -1;
                    comm->state = -1;
                    pthread_mutex_lock(mutex);
                    manager->state = -1;
                    pthread_cond_broadcast(cond);
                    pthread_mutex_unlock(mutex);
                    break;


                // case 2: //cleanup
                //     if (manager->x)
                //     {
                //         delete manager->x;
                //         manager->x = nullptr;
                //     }
                //     manager->state = 0;
                //     break;


                case 1: //recalculate
                    pthread_mutex_lock(mutex);
                    // printf("here\n");
                    manager->nx = comm->nx;
                    manager->ny = comm->ny;
                    manager->f = comm->f;
                    manager->fMax = comm->fMax;
                    manager->outlier = comm->outlier;
                    n = (manager->nx + 1) * (manager->ny + 1);
                    manager->A = new double[MSR_offD_size(manager->nx, manager->ny) +
                        n + 1];
                    if (!(manager->A))
                    {
                        err = err | 1;
                        manager->A = nullptr;
                    }
                    manager->I = new unsigned int[MSR_offD_size(manager->nx, manager->ny) +
                        n + 1];
                    if (!(manager->I))
                    {
                        err = err | 2;
                        manager->I = nullptr;
                    }
                    manager->b = new double[n];
                    if (!(manager->b))
                    {
                        err = err | 4;
                        manager->b = nullptr;
                    }
                    manager->x = new double[n];
                    if (!(manager->x))
                    {
                        err = err | 8;
                        manager->x = nullptr;
                    }
                    manager->r = new double[n];
                    if (!(manager->r))
                    {
                        err = err | 16;
                        manager->r = nullptr;
                    }
                    manager->u = new double[n];
                    if (!(manager->u))
                    {
                        err = err | 32;
                        manager->u = nullptr;
                    }
                    manager->v = new double[n];
                    if (!(manager->v))
                    {
                        err = err | 64;
                        manager->v = nullptr;
                    }
                    manager->buffer = new double[p];
                    if (!(manager->buffer))
                    {
                        err = err | 128;
                        manager->buffer = nullptr;
                    }



                    if (!err)
                    {
                        manager->state = 1;
                        comm->state = 1;
                        pthread_cond_broadcast(cond);
                        pthread_mutex_unlock(mutex);
                        pthread_mutex_unlock(&(comm->guitex));
                        ret = solve(manager, p, t, barrier);
                        // printf("almost\n");
                        pthread_mutex_lock(&(comm->guitex));
                        pthread_mutex_lock(mutex);
                    }


                    if (!(err & 1))
                    {
                        delete[] manager->A;
                    }
                    if (!(err & 2))
                    {
                        delete[] manager->I;
                    }
                    if (!(err & 4))
                    {
                        delete[] manager->b;
                    }
                    if (!(err & 16))
                    {
                        delete[] manager->r;
                    }
                    if (!(err & 32))
                    {
                        delete[] manager->u;
                    }
                    if (!(err & 64))
                    {
                        delete[] manager->v;
                    }
                    if (!(err & 128))
                    {
                        delete[] manager->buffer;
                    }
                    manager->A = nullptr;
                    manager->I = nullptr;
                    manager->b = nullptr;
                    manager->r = nullptr;
                    manager->u = nullptr;
                    manager->v = nullptr;
                    manager->buffer = nullptr;
                    manager->err = err;
                    comm->err = err;
                    comm->iteration = ret;
                    manager->state = 0;
                    comm->state = 0;

                    // comm->nx = manager->nx;
                    // comm->ny = manager->ny;


                    if ((!err) && (ret >= 0))
                    {
                        comm->state = 2;

                        if (comm->x)
                            delete[] comm->x;

                        comm->x = new double[n];
                        memcpy(comm->x, manager->x, n * sizeof(double));
                    }
                    comm->wakeGUI();


                    if (!(err & 8))
                    {
                        delete[] manager->x;
                    }
                    manager->x = nullptr;

                    err = 0;
                    pthread_mutex_unlock(mutex);

                    if (comm->esheodin_flg != 1)
                    {
                        comm->esheodin_flg = 1;
                        pthread_cond_broadcast(&(comm->esheodin));
                    }
                    // printf("where\n");
                    break;


                default:
                    break;
            }

            comm->request = 0;
        }

        delete []comm->x;
        comm->x = nullptr;
        pthread_mutex_unlock(&(comm->guitex));
    }
    else
    {
        pthread_mutex_lock(mutex);
        while (quit != -1)
        {
            manager->count++;
            if (manager->count >= p)
            {
                manager->count = 0;
                manager->ready_flg = 1;
                pthread_cond_broadcast(ready);
                // printf("there %d\n", t);
            }
            pthread_cond_wait(cond, mutex);
            // printf("cccccc %d\n", t);
            switch(manager->state)
            {
                case -1:
                    quit = -1;
                    // printf("fffkkkkk%d\n", t);
                    break;

                case 1:
                    // printf("ddddddd %d\n", t);
                    pthread_mutex_unlock(mutex);
                    // printf("vvvvvvv %d\n", t);
                    solve(manager, p, t, barrier);
                    // printf("almost %d\n", t);
                    pthread_mutex_lock(mutex);

                    
                    // printf("every %d\n", t);

                    break;

                default:
                    break;
            }
        }
        pthread_mutex_unlock(mutex);
    }

    // printf("fff%d\n", t);
    pthread_barrier_wait(barrier);
    if (t == 0)
    {
        pthread_mutex_lock(&(comm->guitex));
        comm->rquit_flg = 1;
        pthread_cond_broadcast(&(comm->rquit));
        // printf("fff\n");
        pthread_mutex_unlock(&(comm->guitex));
    }

    return nullptr;
}






calcThreadManager::calcThreadManager(double _x1, double _x2, double _y1, double _y2,
     int _maxit, double _eps, int _p)
{
    x1 = _x1;
    x2 = _x2;
    y1 = _y1;
    y2 = _y2;
    maxit = _maxit;
    eps = _eps;
    p = _p;
	pthread_barrier_init(&barrier, 0, p);
}


int calcThreadManager::getThreadArgs(thread_args* args, guiThreadComm *comm)
{
    for(int t = 0; t < p; t++)
    {
        args[t].comm = nullptr;
        args[t].manager = this;
        args[t].p = p;
        args[t].t = t;
    }
    args[0].comm = comm;

    return 0;
}


calcThreadManager::~calcThreadManager()
{
    if (A)
        delete[] A;
    if (I)
        delete[] I;
    if (b)
        delete[] b;
    if (x)
        delete[] x;
    if (r)
        delete[] r;
    if (u)
        delete[] u;
    if (v)
        delete[] v;
    if (buffer)
        delete[] buffer;
    A = nullptr;
    I = nullptr;
    b = nullptr;
    x = nullptr;
    r = nullptr;
    u = nullptr;
    v = nullptr;
    buffer = nullptr;
	pthread_barrier_destroy(&barrier);
}






guiThreadComm::guiThreadComm() : wakeEventType((QEvent::Type)QEvent::registerEventType())
{

}


int guiThreadComm::wakeGUI()
{
    QEvent *event = new QEvent(wakeEventType);
    if (app)
    {
        // printf("I send\n");
        QCoreApplication::postEvent(app, event);
        // printf("I go\n");
        return 0;
    }
    return -1;
}


int guiThreadComm::attachApp(QObject *_app)
{
    app = _app;
    return 0;
}


guiThreadComm::~guiThreadComm()
{
    if (x)
        delete[] x;
    x = nullptr;
}