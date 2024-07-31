#ifndef __WINDOW_H__
#define __WINDOW_H__

#include <QtWidgets/QtWidgets>
#include <QOpenGLWidget>
#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QSurfaceFormat>
#include <QColor>

#include <iostream>
#include <stdio.h>

#include "threadManager.h"

#define W_RES_DEFAULT 128
#define B_RES_DEFAULT 128

struct v3
{
    double x;
    double y;
    double z;
};


class approxWindow : public QOpenGLWidget, protected QOpenGLFunctions
{
    // Q_OBJECT

    private:
        class renderSurface
        {
            public:
                approxWindow *owner;
                double *vertexArray;
                unsigned int *indexArray;
                unsigned int wRes = W_RES_DEFAULT;
                unsigned int bRes = B_RES_DEFAULT;
                unsigned int err = 0;
                double zMin = 0., zMax = 0.;

            public:
                renderSurface(approxWindow* _owner);
                unsigned int update();
                unsigned int vertexCount();
                unsigned int indexCount();
                int update_fBounds(double &fMin, double &fMax);
                int update_fBounds(double &fMin, double &fMax, double (*f)(double, double));
                ~renderSurface();

        };


    private:
        unsigned int xPt, yPt;
        int s, p, r;
        int iterations = 0;
        double x1, x2, y1, y2;
        int funcID, graphID;
        const char *funcName;
        double (*f)(double, double);
        double fMin, fMax;
        double graphMin = 0., graphMax = 0.;
        // float rightPlane = 0.f, leftPlane = 0.f, topPlane = 0.f, bottomPlane = 0.f;
        // float nearPlane = 0.f, farPlane = 0.f;
        QShortcut* hotkeys[10];
        approxWindow::renderSurface *surface;

        int isWaiting = 0;
        int rescale = 0;

        double* X = nullptr;

        float axisVertexArray[12] = 
            {0.f, 0.f, 0.f,
             1.f, 0.f, 0.f,
             0.f, 0.f, 1.f,
             0.f, 1.f, 0.f};
        unsigned int axisIndexArray[6] = {0, 1, 0, 2, 0, 3};
        unsigned int ssIndexArray[6] = {0, 1, W_RES_DEFAULT, 1, W_RES_DEFAULT, W_RES_DEFAULT + 1};

    private:
        QOpenGLVertexArrayObject surfaceVAO;
        GLuint surfaceVBO;
        GLuint surfaceEBO;

        QOpenGLVertexArrayObject axisVAO;
        GLuint axisVBO;
        GLuint axisEBO;

        int translateSurfaceMatLoc = 0;
        int translateAxisMatLoc = 0;
        int rotateSurfaceMatLoc = 0;
        int rotateAxisMatLoc = 0;
        int scaleSurfaceMatLoc = 0;
        int scaleAxisMatLoc = 0;
        int projectionSurfaceMatLoc = 0;
        int projectionAxisMatLoc = 0;
        int minLoc = 0;
        int maxLoc = 0;
        QMatrix4x4 translate;
        QMatrix4x4 rotate;
        QMatrix4x4 scale;
        QMatrix4x4 scaleAxis;
        QMatrix4x4 projection;

        QOpenGLShaderProgram *surfaceShaderProg;
        QOpenGLShaderProgram *axisShaderProg;

    private: 
        guiThreadComm *comm;

    public:
        approxWindow(QWidget *parent, guiThreadComm *_comm, int /*k*/,
     double _x1, double _x2, double _y1, double _y2, int nx, int ny);
        void* start();
        void bboxScale();
        QSize minimumSizeHint() const;
        QSize sizeHint() const;
        double approxFunc(double x, double y);
        double errorFunc(double x, double y);
        double renderFunc(double x, double y);
        void customEvent(QEvent *event);
        int requestApproxUpdate(unsigned int nx, unsigned int ny);
        int surfaceUpdate();
        void getInfoString();
        QString getInfoStringScreen();
        ~approxWindow();


    public slots:
        void func_cycle();
        void graph_cycle();
        void s_inc();
        void s_decr();
        void n_inc();
        void n_decr();
        void p_inc();
        void p_decr();
        void r_incr();
        void r_decr();


    protected:
        void initializeGL();
        void resizeGL(int width, int height);
        void paintGL();




};

#endif