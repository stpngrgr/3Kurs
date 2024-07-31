#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QShortcut>
#include <QKeySequence>

#include "window.h"

int main(int argc, char **argv)
{
    int err;
    QShortcut* quit;

    QApplication app(argc, argv);

    QMainWindow *window = new QMainWindow;
    ApproxWindow *aprxwndw = new ApproxWindow(window);
    err = aprxwndw->parse_command_line(argc, argv);

    if (err > 0)
        printf("No command line arguments passed, default arguments used instead\nUsage : %s a b n k\n", argv[0]);
    if (err < 0)
        switch (-err)
        {
            case 1:
                printf("Bad segment\nUsage : %s a b n k\n", argv[0]);
                return -1;
                break;
            case 2:
                printf("Bad n\nUsage : %s a b n k\n", argv[0]);
                return -1;
                break;
            case 3:
                printf("Bad k\nUsage : %s a b n k\n", argv[0]);
                return -1;
                break;
            case 5:
                printf("Wrong argument count\nUsage : %s a b n k\n", argv[0]);
                return -1;
                break;
            default:
                printf("Unknown initialization error\nUsage : %s a b n k\n", argv[0]);
                return -1;
                break;
        }

    window->setCentralWidget(aprxwndw);
    window->setWindowTitle ("Test");

    quit = new QShortcut(QKeySequence("q"), window, SLOT(close()));
    // printf("quit : %p\n", (void*)quit);

    window->show();
    app.exec();
    // printf("quit : %p\n", (void*)quit);
    if (quit)
    {
        delete quit;
        quit = nullptr;
    }
    // delete aprxwndw;
    if (window)
    {
        delete window;
        window = nullptr;
    }
    return 0;
}