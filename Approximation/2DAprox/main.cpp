#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#define __MAXIT__ 10000
#define __EPS__ 1.e-14


#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QShortcut>
#include <QKeySequence>


#include "window.h"
#include "threadManager.h"






int readArea(char *filename, double *data)
{
    FILE *f;
    char c;
    int err;
    double temp;
    int dataCount = 0;

    f = fopen(filename, "r");
    if (!f)
        return -1;

    while (!feof(f))
    {
        while (fscanf(f, "%lf", &temp) == 1)
        {
            if (dataCount < 4)
                data[dataCount++] = temp;
            else
            {
                fclose(f);
                return -2;
            }
        }

        if (!feof(f)) 
        {
            err = fscanf(f, "%c", &c);
            // printf("%c", c);
            if ((err == 1) && (c == '#'))
            {
                while ((fscanf(f, "%c", &c) == 1) && (c != '\n')) 
                { }
            }
        }
    }

    if (dataCount != 4)
    {
        fclose(f);
        return -3;
    }

    fclose(f);
    return 0;
}






int main(int argc, char **argv)
{
    int err = 0;
    double data[4];
    int maxit = __MAXIT__;
    double eps = __EPS__;
    int p;
    int nx, ny;
    int k;
	pthread_t* threadIds;
    thread_args* threadArgs;

    QShortcut* quit;

    QApplication app(argc, argv);

    //process input
    if (argc != 7)
    {
        printf("Usage : %s filename nx ny k eps p\n", argv[0]);
        return -1;
    }
    if ((sscanf(argv[2], "%d", &nx) != 1) || (nx < 0))
    {
        printf("Bad nx\n");
        printf("Usage : %s filename nx ny k eps p\n", argv[0]);
        return -1;
    }
    if ((sscanf(argv[3], "%d", &ny) != 1) || (ny < 0))
    {
        printf("Bad ny\n");
        printf("Usage : %s filename nx ny k eps p\n", argv[0]);
        return -1;
    }
    if ((sscanf(argv[4], "%d", &k) != 1) || (k < 0) || (k > 7))
    {
        printf("Bad k\n");
        printf("Usage : %s filename nx ny k eps p\n", argv[0]);
        return -1;
    }
    if ((sscanf(argv[5], "%lf", &eps) != 1) || (eps < 1.e-15))
    {
        printf("Bad eps\n");
        printf("Usage : %s filename nx ny k eps p\n", argv[0]);
        return -1;
    }
    if ((sscanf(argv[6], "%d", &p) != 1) || (p < 0))
    {
        printf("Bad p\n");
        printf("Usage : %s filename nx ny k eps p\n", argv[0]);
        return -1;
    }

    err = readArea(argv[1], data);
    if (err)
    {
        printf("File reading error\n");
        return -1;
    }


    guiThreadComm comm;
    calcThreadManager manager(data[0], data[2], data[3], data[1], maxit, eps, p);

    QMainWindow *window = new QMainWindow;
    approxWindow *aprxwndw = new approxWindow(window, &comm, k, data[0], data[2], data[3], data[1], nx, ny);

    comm.attachApp(aprxwndw);


    threadIds = new pthread_t[p];
    if (threadIds == nullptr)
    {
        printf("threadIds creation failed\n");
        delete aprxwndw;
        delete window;
        return -1;
    }
    threadArgs = new thread_args[p];
    if (threadIds == nullptr)
    {
        printf("threadArgs creation failed\n");
        delete aprxwndw;
        delete window;
        delete[] threadIds;
        return -1;
    }
    manager.getThreadArgs(threadArgs, &comm);
    for (int i = 0; i < p; i++)
    {
        err = pthread_create(threadIds, 0, thread_do, (void *)(threadArgs + i));
        if (err)
        {
            printf(" pthread_create failed%d\n", i);
            delete aprxwndw;
            delete window;
            delete[] threadIds;
            delete[] threadArgs;
            return -1;
        }
    }



    window->setCentralWidget(aprxwndw);
    window->setWindowTitle ("Test");

    quit = new QShortcut(QKeySequence("q"), window, SLOT(close()));

    aprxwndw->start();

    // printf("doged it\n");

    // printf("g\n");

    window->show();
    app.exec();

    pthread_mutex_lock(&(comm.guitex));
    while (!comm.ready_flg)
        pthread_cond_wait(&(comm.ready), &(comm.guitex));

    comm.request = -1;
    comm.rquit_flg = 0;
    pthread_cond_broadcast(&(comm.wakeupcall));
    // printf("fffddddddddddf\n");

    while (!comm.rquit_flg)
        pthread_cond_wait(&(comm.rquit), &(comm.guitex));
    
    pthread_mutex_unlock(&(comm.guitex));


    delete quit;
    delete aprxwndw;
    delete window;
    delete[] threadIds;
    delete[] threadArgs;
    // printf("exit : gui\n");
    return 0;
}