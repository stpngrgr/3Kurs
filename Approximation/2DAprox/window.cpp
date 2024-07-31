#include <math.h>
#include "window.h"

#define pi 3.1415






// static void display(double *x, unsigned int n)
// {
//     printf("\n\n ---------------- \n\n");
//     for (unsigned int i = 0; i < n; i++)
//     {
//         printf("%e\n", x[i]);
//     }
//     printf("\n\n ---------------- \n\n");
// }







static const char *surfaceVSSource = 
    "#version 330 core\n"
    "layout (location = 0)in vec3 vertex;\n"
    "uniform mat4 rotate;\n"
    "uniform mat4 scale;\n"
    "uniform mat4 translate;\n"
    "uniform mat4 projection;\n"
    "out float h;\n"
    "void main() {\n"
    // "   gl_Position = vec4(vertex.x, vertex.y, vertex.z, 1.0);\n"
    "   vec4 v = projection * translate * rotate * scale * vec4(vertex.x, vertex.y, vertex.z, 1.0);\n"
    "   gl_Position = v;\n"
    "   h = vertex.y;\n"
    "}\n";

static const char *axisVSSource = 
    "#version 330 core\n"
    "layout (location = 0)in vec3 vertex;\n"
    "uniform mat4 rotate;\n"
    "uniform mat4 scale;\n"
    "uniform mat4 translate;\n"
    "uniform mat4 projection;\n"
    "void main() {\n"
    // "   gl_Position = vec4(vertex.x, vertex.y, vertex.z, 1.0);\n"
    "   gl_Position = projection * translate * rotate * scale * vec4(vertex.x, vertex.y, vertex.z, 1.0);\n"
    "}\n";

static const char *surfaceFSSource = 
    "#version 330 core\n"
    "in float h;"
    "out vec4 FragColour;"
    "uniform float min;\n"
    "uniform float max;\n"
    "void main() {\n"
    "   float c = (h - min) / (max - min);\n"
    "   FragColour = vec4(c * 1.3, 0.2, 0.8, 1.0);\n"
    "}\n";

static const char *axisFSSource = 
    "#version 330 core\n"
    "out vec4 FragColour;"
    "void main() {\n"
    "   FragColour = vec4(0.0, 0.0, 0.0, 0.0);\n"
    "}\n";




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







double f0(double /*x*/, double /*y*/) //f = 1;
{
    return 1.;
}

double f1(double x, double /*y*/) //f = x;
{
    return x;
}

double f2(double /*x*/, double y) //f = y;
{
    return y;
}

double f3(double x, double y) //f = x + y;
{
    return x + y;
}

double f4(double x, double y) //f = sqrt(x * x + y * y);
{
    return sqrt(x * x + y * y);
}

double f5(double x, double y) //f = x * x + y * y;
{
    return x * x + y * y;
}

double f6(double x, double y) //f = exp(x * x - y * y);
{
    return exp(x * x - y * y);
}

double f7(double x, double y) //f = 1. / (25.(x * x + y * y) + 1.);
{
    return 1. / (25. * (x * x + y * y) + 1.);
}




double approxWindow::approxFunc(double x, double y)
{
    double x_0 = x1, y_0 = y1, z_0;
    double x_1 = x2, y_1 = y2, z_1;
    double x_2, y_2, z_2;

    double delta_x = (x2 - x1) / (xPt - 1);
    double delta_y = (y2 - y1) / (yPt - 1);

    unsigned int il = 0, ir = xPt - 1, jl = 0, jr = yPt - 1;
    unsigned int i, j, l;

    double a[9];
    double temp;

    while (ir > il + 1)
    {
        i = (ir + il) / 2;
        z_2 = x1 + i * delta_x;
        if (z_2 < x)
        {
            x_0 = z_2;
            il = i;
        }
        else
        {
            x_1 = z_2;
            ir = i;
        }
    }

    while (jr > jl + 1)
    {
        j = (jr + jl) / 2;
        z_2 = y1 + j * delta_y;
        if (z_2 < y)
        {
            y_0 = z_2;
            jl = j;
        }
        else
        {
            y_1 = z_2;
            jr = j;
        }
    }

    ij2l(xPt - 1, yPt - 1, il, jl, l);
    z_0 = X[l];
    ij2l(xPt - 1, yPt - 1, ir, jr, l);
    z_1 = X[l];
    if ((x - x_0) > (y - y_0))
    {
        i = ir;
        j = jl;
    }
    else
    {
        i = il;
        j = jr;
    }

    x_2 = x1 + i * delta_x;
    y_2 = y1 + j * delta_y;
    ij2l(xPt - 1, yPt - 1, i, j, l);
    z_2 = X[l];

    a[0] = x - x_0;
    a[1] = y - y_0;
    a[2] = 0.;
    a[3] = x_1 - x_0;
    a[4] = y_1 - y_0;
    a[5] = z_1 - z_0;
    a[6] = x_2 - x_0;
    a[7] = y_2 - y_0;
    a[8] = z_2 - z_0;

    a[2] -= a[0] * (a[4] * a[8] - a[5] * a[7]);
    a[2] += a[1] * (a[3] * a[8] - a[5] * a[6]);
    temp = (a[3] * a[7] - a[4] * a[6]);
    a[2] += z_0 * temp;

    // printf("(%e, %e) : { (%e, %e, %e) , (%e, %e, %e) , (%e, %e, %e) }\n",
    //     x, y, x_0, y_0, z_0, x_1, y_1, z_1, x_2, y_2, z_2);
    // printf("{ (%d, %d) , (%d, %d) , (%d, %d) } : %e\n",
    //     il, jl, ir, jr, i, j, a[2] / temp);
    // printf ("-  %e  +  %e  +  %e\n", a[0] * (a[4] * a[8] - a[5] * a[7]),
    //      a[1] * (a[3] * a[8] - a[5] * a[6]), z_0 * temp);

    return a[2] / temp;
}





double approxWindow::errorFunc(double x, double y)
{
    return fabs(f(x, y) - approxFunc(x, y));
}







void approxWindow::func_cycle()
{
    if (!isWaiting)
    {
        funcID = (funcID + 1) % 8;

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
                funcName = "k = 2 | f(x) = y";
                break;
            case 3:
                f = f3;
                funcName = "k = 3 | x + y";
                break;
            case 4:
                f = f4;
                funcName = "k = 4 | f(x) = sqrt(x^2 + y^2)";
                break;
            case 5:
                f = f5;
                funcName = "k = 5 | f(x) = x^2 + y^2";
                break;
            case 6:
                f = f6;
                funcName = "k = 6 | f(x) = e^(x^2 - y^2)";
                break;
            case 7:
                f = f7;
                funcName = "k = 6 | f(x) = 1/(25(x^2 + y^2) + 1)";
                break;
            default:
                f = f0;
                funcName = "k = 0 | f(x) = 1";
                perror("func_cycle default case executed\n");
                break;
        }

        rescale = 1;
        requestApproxUpdate(xPt, yPt);
    }
}


void approxWindow::graph_cycle()
{
    graphID = (graphID + 1) % 3;
        
    if (graphID != 0)
        surface->update_fBounds(graphMin, graphMax);
    else
    {
        graphMin = fMin;
        graphMax = fMax;
    }
    surfaceUpdate();
}


void approxWindow::s_inc()
{
    s++;
    bboxScale();
    update();
}


void approxWindow::s_decr()
{
    s--;
    if (s < 0)
        s = 0;
    bboxScale();
    update();
}


void approxWindow::n_inc()
{
    unsigned int _nx = xPt, _ny = yPt;
    if (!isWaiting)
    {
        _nx *= 2;
        _ny *= 2;
        requestApproxUpdate(_nx, _ny);
    }
}


void approxWindow::n_decr()
{
    unsigned int _nx = xPt, _ny = yPt;
    if (!isWaiting)
    {
        // printf("\n\n ndecr\n\n");
        _nx /= 2;
        _ny /= 2;
        if ((_nx >= 2) && (_ny >= 2))
        {
            requestApproxUpdate(_nx, _ny);
        }
    }
}


void approxWindow::p_inc()
{
    if (!isWaiting)
    {
        p++;
        requestApproxUpdate(xPt, yPt);
    }
}


void approxWindow::p_decr()
{
    if (!isWaiting)
    {
        p--;
        requestApproxUpdate(xPt, yPt);
    }
}

void approxWindow::r_incr()
{
    r = (r + 1) % 24;
    getInfoString();
    update();
}


void approxWindow::r_decr()
{
    r -= 1;
    if (r < 0)
        r = 23;
    getInfoString();
    update();
}





approxWindow::approxWindow(QWidget *parent, guiThreadComm *_comm, int k,
     double _x1, double _x2, double _y1, double _y2, int nx, int ny) : QOpenGLWidget (parent)
{
    p = 0;
    s = 0;
    r = 0;
    x1 = _x1;
    x2 = _x2;
    // y1 = 5.;
    // y2 = -5.;
    y1 = _y1;
    y2 = _y2;
    xPt = (unsigned int)nx;
    yPt = (unsigned int)ny;
    funcID = k;
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
                funcName = "k = 2 | f(x) = y";
                break;
            case 3:
                f = f3;
                funcName = "k = 3 | x + y";
                break;
            case 4:
                f = f4;
                funcName = "k = 4 | f(x) = sqrt(x^2 + y^2)";
                break;
            case 5:
                f = f5;
                funcName = "k = 5 | f(x) = x^2 + y^2";
                break;
            case 6:
                f = f6;
                funcName = "k = 6 | f(x) = e^(x^2 - y^2)";
                break;
            case 7:
                f = f7;
                funcName = "k = 6 | f(x) = 1/(25(x^2 + y^2) + 1)";
                break;
            default:
                f = f0;
                funcName = "k = 0 | f(x) = 1";
                perror("func_cycle default case executed\n");
                break;
        }
    graphID = 0;
    comm = _comm;
    X = new double[xPt * yPt];


    hotkeys[0] = new QShortcut(QKeySequence("0"), (approxWindow*)this);
    connect(hotkeys[0], &QShortcut::activated, this, &approxWindow::func_cycle);
    hotkeys[1] = new QShortcut(QKeySequence("1"), (approxWindow*)this);
    connect(hotkeys[1], &QShortcut::activated, this, &approxWindow::graph_cycle);
    hotkeys[2] = new QShortcut(QKeySequence("2"), (approxWindow*)this);
    connect(hotkeys[2], &QShortcut::activated, this, &approxWindow::s_inc);
    hotkeys[3] = new QShortcut(QKeySequence("3"), (approxWindow*)this);
    connect(hotkeys[3], &QShortcut::activated, this, &approxWindow::s_decr);
    hotkeys[4] = new QShortcut(QKeySequence("4"), (approxWindow*)this);
    connect(hotkeys[4], &QShortcut::activated, this, &approxWindow::n_inc);
    hotkeys[5] = new QShortcut(QKeySequence("5"), (approxWindow*)this);
    connect(hotkeys[5], &QShortcut::activated, this, &approxWindow::n_decr);
    hotkeys[6] = new QShortcut(QKeySequence("6"), (approxWindow*)this);
    connect(hotkeys[6], &QShortcut::activated, this, &approxWindow::p_inc);
    hotkeys[7] = new QShortcut(QKeySequence("7"), (approxWindow*)this);
    connect(hotkeys[7], &QShortcut::activated, this, &approxWindow::p_decr);
    hotkeys[8] = new QShortcut(QKeySequence("8"), (approxWindow*)this);
    connect(hotkeys[8], &QShortcut::activated, this, &approxWindow::r_incr);
    hotkeys[9] = new QShortcut(QKeySequence("9"), (approxWindow*)this);
    connect(hotkeys[9], &QShortcut::activated, this, &approxWindow::r_decr);


    surface = new renderSurface(this);


    translate.setToIdentity();
    rotate.setToIdentity();
    rotate.lookAt(QVector3D(1.f, 1.f, 1.f), QVector3D(0.f, 0.f, 0.f), QVector3D(0.f, 1.f, 0.f));
    surface->update_fBounds(graphMin, graphMax);
    bboxScale();
}




void* approxWindow::start()
{
    pthread_mutex_lock(&(comm->guitex));
    // printf("me too!\n");
    // comm->count++;
    // if (comm->count >= 2)
    // {
    //     comm->count = 0;
    //     comm->rquit_flg = 1;
    //     pthread_cond_broadcast(&(comm->rquit));
    //     printf("see what i can do!\n");
    // }
    // else
    // {
    //     while (!comm->ready_flg)
    //         pthread_cond_wait(&(comm->ready), &(comm->guitex));
    // }
    // printf("wonderful!\n");
    while (!comm->ready_flg)
        pthread_cond_wait(&(comm->ready), &(comm->guitex));
    comm->ready_flg = 0;
    pthread_mutex_unlock(&(comm->guitex));
    // printf("a\n");
    requestApproxUpdate(xPt, yPt);
    pthread_mutex_lock(&(comm->guitex));
    while (!comm->esheodin_flg)
        pthread_cond_wait(&(comm->esheodin), &(comm->guitex));
    pthread_mutex_unlock(&(comm->guitex));
    return nullptr;
}
        




int approxWindow::requestApproxUpdate(unsigned int nx, unsigned int ny)
{
    // int upd = 0;

    // printf("brother!\n");
    pthread_mutex_lock(&(comm->guitex));
    // printf("is it real? %d\n", comm->state);
    // if (comm->state == 2)
    // {
    //     xPt = comm->nx + 1;
    //     yPt = comm->ny + 1;
    //     if (X)
    //         delete[] X;
    //     X = new double[xPt * yPt];
    //     memcpy(X, comm->x, xPt * yPt * sizeof(double));
    //     comm->state = 0;
    //     upd = 1;
    // }
    if (comm->state == 0)
    {
        // printf("or is it not?\n");
        comm->nx = nx - 1;
        comm->ny = ny - 1;
        comm->f = f;
        comm->outlier = p;

        surface->update_fBounds(fMin, fMax, f);

        comm->fMax = (fabs(fMax) > fabs(fMin) ? fabs(fMax) : fabs(fMin));
        comm->request = 1;
        pthread_cond_broadcast(&(comm->wakeupcall));
        isWaiting = 1;
        // printf("isnt it?\n");
    }
    // if (comm->state == -1)
    // {
    //     q = 1;
    // }
    pthread_mutex_unlock(&(comm->guitex));

    // if (upd)
    //     surfaceUpdate();

    // if (q)
    //     parent->

    return 0;
}





void approxWindow::bboxScale()
{
    float maxScale = 0.f, temp = 0.f;
    double z = (fabs(graphMax) > fabs(graphMin) ? fabs(graphMax) : fabs(graphMin));

    temp = (float)(sqrt(x1 * x1 + y1 * y1) + z);
    if (temp > maxScale)
    {
        maxScale = temp;
    }
    temp = (float)(sqrt(x1 * x1 + y2 * y2) + z);
    if (temp > maxScale)
    {
        maxScale = temp;
    }
    temp = (float)(sqrt(x2 * x2 + y1 * y1) + z);
    if (temp > maxScale)
    {
        maxScale = temp;
    }
    temp = (float)(sqrt(x2 * x2 + y2 * y2) + z);
    if (temp > maxScale)
    {
        maxScale = temp;
    }

    maxScale /= (float)pow(2., s);

    scale.setToIdentity();
    scaleAxis.setToIdentity();
    scaleAxis.scale(maxScale, maxScale, maxScale);
    projection.setToIdentity();  
    // projection.ortho(-10.f, 10.f, -10.f, 10.f, -10.f, 10.f);  
    projection.ortho(-0.8 * maxScale, 0.8 * maxScale, -0.8 * maxScale, 0.8 * maxScale, -10., 10.);

    return;
}






QSize approxWindow::minimumSizeHint() const
{
    return QSize(100, 100);
}


QSize approxWindow::sizeHint() const
{
    return QSize(1000, 1000);
}








double approxWindow::renderFunc(double x, double y)
{
    switch(graphID)
    {
        case 0 :
            return f(x, y);
        case 1 :
            return approxFunc(x, y);
        case 2 :
            return errorFunc(x, y);
        default :
            printf("Default case\n");
            return  f(x, y);
    }

    return  -1.;
}





void approxWindow::initializeGL()
{
    bool t1 = false, t2 = false;

    makeCurrent();
    initializeOpenGLFunctions();
    glViewport(0, 0, width(), height());
    glEnable(GL_DEPTH_TEST);
    glClearColor(1., 1., 1., 1.);

    // qDebug() << (const char*)(QOpenGLContext::currentContext()->functions()->glGetString(GL_VERSION));
    // qDebug() << (const char*)(QOpenGLContext::currentContext()->functions()->glGetString(GL_SHADING_LANGUAGE_VERSION));


    rotate.setToIdentity();
    rotate.rotate(r * 15.f, 0.f, 1.f, 0.f);
    rotate.lookAt(QVector3D(1.f, 1.f, 1.f), QVector3D(0.f, 0.f, 0.f), QVector3D(0.f, 1.f, 0.f));


    surfaceShaderProg = new QOpenGLShaderProgram;
    t1 = surfaceShaderProg->addShaderFromSourceCode(QOpenGLShader::Vertex, surfaceVSSource);
    t2 = surfaceShaderProg->addShaderFromSourceCode(QOpenGLShader::Fragment, surfaceFSSource);
    if ((!t1) || (!t2))
        printf("Shader linkage failed\n");
    surfaceShaderProg->bindAttributeLocation("vertex", 0);
    surfaceShaderProg->link();
    surfaceShaderProg->bind();
    translateSurfaceMatLoc = surfaceShaderProg->uniformLocation("translate");
    rotateSurfaceMatLoc = surfaceShaderProg->uniformLocation("rotate");
    scaleSurfaceMatLoc = surfaceShaderProg->uniformLocation("scale");
    projectionSurfaceMatLoc = surfaceShaderProg->uniformLocation("projection");
    minLoc = surfaceShaderProg->uniformLocation("min");
    maxLoc = surfaceShaderProg->uniformLocation("max");

    surfaceVAO.create();
    // QOpenGLVertexArrayObject::Binder vaoBinder(&surfaceVAO);
    surfaceVAO.bind();

    glGenBuffers(1, &surfaceVBO);
    glBindBuffer(GL_ARRAY_BUFFER, surfaceVBO);
    glBufferData(GL_ARRAY_BUFFER, surface->vertexCount() * sizeof(double), surface->vertexArray, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), 0);

    glGenBuffers(1, &surfaceEBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, surfaceEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, surface->indexCount() * sizeof(unsigned int), surface->indexArray, GL_STATIC_DRAW);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * sizeof(unsigned int), ssIndexArray, GL_STATIC_DRAW);

    surfaceVAO.release();
    surfaceShaderProg->release();


    axisShaderProg = new QOpenGLShaderProgram;
    t1 = axisShaderProg->addShaderFromSourceCode(QOpenGLShader::Vertex, axisVSSource);
    t2 = axisShaderProg->addShaderFromSourceCode(QOpenGLShader::Fragment, axisFSSource);
    if ((!t1) || (!t2))
        printf("Shader linkage failed\n");
    axisShaderProg->bindAttributeLocation("vertex", 0);
    axisShaderProg->link();
    axisShaderProg->bind();
    translateAxisMatLoc = axisShaderProg->uniformLocation("translate");
    rotateAxisMatLoc = axisShaderProg->uniformLocation("rotate");
    scaleAxisMatLoc = axisShaderProg->uniformLocation("scale");
    projectionAxisMatLoc = axisShaderProg->uniformLocation("projection");

    axisVAO.create();
    // QOpenGLVertexArrayObject::Binder vaoBinder(&axisVAO);
    axisVAO.bind();

    glGenBuffers(1, &axisVBO);
    glBindBuffer(GL_ARRAY_BUFFER, axisVBO);
    glBufferData(GL_ARRAY_BUFFER, 12 * sizeof(float), axisVertexArray, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), 0);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &axisEBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, axisEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * sizeof(unsigned int), axisIndexArray, GL_STATIC_DRAW);

    axisVAO.release();
    axisShaderProg->release();

}


void approxWindow::resizeGL(int width, int height)
{
    glViewport(0, 0, width, height);
}


void approxWindow::paintGL()
{
    QPainter painter(this);
    QPen black(Qt::black, 0, Qt::SolidLine);

    painter.setPen(black);
    painter.drawText (100, 100, getInfoStringScreen());

    painter.beginNativePainting();
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // glEnable(GL_CULL_FACE);
    glViewport(0, 0, width(), height());

    rotate.setToIdentity();
    rotate.lookAt(QVector3D(1.f, 1.f, 1.f), QVector3D(0.f, 0.f, 0.f), QVector3D(0.f, 1.f, 0.f));
    // rotate.rotate(r * 15, 0.f, 1.f, 0.f);
    rotate.rotate((r % 24) * 15, 0.f, 1.f, 0.f);

    axisShaderProg->bind();
    axisShaderProg->setUniformValue(translateAxisMatLoc, translate);
    axisShaderProg->setUniformValue(rotateAxisMatLoc, rotate);
    axisShaderProg->setUniformValue(scaleAxisMatLoc, scaleAxis);
    axisShaderProg->setUniformValue(projectionAxisMatLoc, projection);
    axisVAO.bind();
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, axisEBO);

    // glEnableVertexAttribArray(0);
    glDrawElements(GL_LINES, 6, GL_UNSIGNED_INT, 0);
    // std::cout << std::hex << glGetError() << std::endl;

    axisVAO.release();
    axisShaderProg->release();


    surfaceShaderProg->bind();
    surfaceShaderProg->setUniformValue(translateSurfaceMatLoc, translate);
    surfaceShaderProg->setUniformValue(rotateSurfaceMatLoc, rotate);
    surfaceShaderProg->setUniformValue(scaleSurfaceMatLoc, scale);
    surfaceShaderProg->setUniformValue(projectionSurfaceMatLoc, projection);
    surfaceShaderProg->setUniformValue(minLoc, (float)(graphMin));
    surfaceShaderProg->setUniformValue(maxLoc, (float)(graphMax));
    surfaceVAO.bind();
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, surfaceEBO);

    // glEnableVertexAttribArray(0);
    // printf("coca\n");
    glDrawElements(GL_TRIANGLES, surface->indexCount(), GL_UNSIGNED_INT, 0);
    // glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    // std::cout << std::hex << glGetError() << std::endl;

    surfaceVAO.release();
    surfaceShaderProg->release();

    painter.endNativePainting();

}




void approxWindow::customEvent(QEvent *event)
{
    int upd = 0;

    if (event->type() == comm->wakeEventType)
    {
        pthread_mutex_lock(&(comm->guitex));
        if (comm->state == 2)
        {
            // printf("cocka\n");
            upd = 1;
            xPt = comm->nx + 1;
            yPt = comm->ny + 1;
            if (X)
                delete[] X;
            X = new double[xPt * yPt];
            memcpy(X, comm->x, xPt * yPt * sizeof(double));
            comm->state = 0;
            iterations = comm->iteration;
            if (graphID != 0)
                surface->update_fBounds(graphMin, graphMax);
            else
            {
                graphMin = fMin;
                graphMax = fMax;
            }
            // display(X, xPt * yPt);
        }
        else if (comm->err)
        {
            printf("error = %d\n", comm->err);
            comm->state = 0;
        }
        else if (comm->iteration < 0)
        {
            printf("Max iterations exceeded\n");
        }
        isWaiting = 0;
        pthread_mutex_unlock(&(comm->guitex));

        if (upd)
            surfaceUpdate();

        return;
    }

    // printf("what???\n");
    return;
}




int approxWindow::surfaceUpdate()
{
    surface->update();

    surfaceVAO.bind();

    glBindBuffer(GL_ARRAY_BUFFER, surfaceVBO);
    glBufferData(GL_ARRAY_BUFFER, surface->vertexCount() * sizeof(double), surface->vertexArray, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 3 * sizeof(double), 0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, surfaceEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, surface->indexCount() * sizeof(unsigned int), surface->indexArray, GL_STATIC_DRAW);

    surfaceVAO.release();

    bboxScale();

    getInfoString();
    update();
    return 0;
}





void approxWindow::getInfoString()
{
    printf("\n\n ---------------------- \n\n");
    printf("%s | graph = %d\n", funcName, graphID);
    printf("nx = %u | ny = %u\n", xPt, yPt);
    printf("[F_m, F_M] = [%e, %e]\n", graphMin, graphMax);
    printf("s = %d | p = %d\n", s, p);
    printf("rotation = %d\n", (r % 24) * 15);
    printf("iterations = %d\n", iterations);
    printf("\n\n ---------------------- \n\n");
    // display(X, xPt * yPt);
}


QString approxWindow::getInfoStringScreen()
{
    QString infoStr = QString("%1 | GraphID = %2 | nx = %3, ny = %4 | [F_m, F_M] = [%5, %6] | s = %7 | p = %8 | r = %9").arg(
        funcName, QString::number(graphID), QString::number(xPt), QString::number(yPt), 
        QString::number(graphMin), QString::number(graphMax), QString::number(s), 
        QString::number(p), QString::number( (r % 24) * 15));

    return infoStr;
}




approxWindow::~approxWindow()
{
    axisVAO.destroy();
    glDeleteBuffers(1, &axisVBO);
    glDeleteBuffers(1, &axisEBO);
    delete axisShaderProg;
    surfaceVAO.destroy();
    glDeleteBuffers(1, &surfaceVBO);
    glDeleteBuffers(1, &surfaceEBO);
    delete surfaceShaderProg;
    doneCurrent();
    delete surface;
    for (int i = 0; i < 10; i++)
    {
        delete hotkeys[i];
    }
    if (X)
        delete[] X;
    X = nullptr;
}