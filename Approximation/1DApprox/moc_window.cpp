/****************************************************************************
** Meta object code from reading C++ file 'window.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "window.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'window.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_ApproxWindow_t {
    QByteArrayData data[11];
    char stringdata0[82];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ApproxWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ApproxWindow_t qt_meta_stringdata_ApproxWindow = {
    {
QT_MOC_LITERAL(0, 0, 12), // "ApproxWindow"
QT_MOC_LITERAL(1, 13, 10), // "func_cycle"
QT_MOC_LITERAL(2, 24, 0), // ""
QT_MOC_LITERAL(3, 25, 11), // "graph_cycle"
QT_MOC_LITERAL(4, 37, 5), // "s_inc"
QT_MOC_LITERAL(5, 43, 6), // "s_decr"
QT_MOC_LITERAL(6, 50, 5), // "n_inc"
QT_MOC_LITERAL(7, 56, 6), // "n_decr"
QT_MOC_LITERAL(8, 63, 5), // "p_inc"
QT_MOC_LITERAL(9, 69, 6), // "p_decr"
QT_MOC_LITERAL(10, 76, 5) // "close"

    },
    "ApproxWindow\0func_cycle\0\0graph_cycle\0"
    "s_inc\0s_decr\0n_inc\0n_decr\0p_inc\0p_decr\0"
    "close"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ApproxWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   59,    2, 0x0a /* Public */,
       3,    0,   60,    2, 0x0a /* Public */,
       4,    0,   61,    2, 0x0a /* Public */,
       5,    0,   62,    2, 0x0a /* Public */,
       6,    0,   63,    2, 0x0a /* Public */,
       7,    0,   64,    2, 0x0a /* Public */,
       8,    0,   65,    2, 0x0a /* Public */,
       9,    0,   66,    2, 0x0a /* Public */,
      10,    0,   67,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void ApproxWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ApproxWindow *_t = static_cast<ApproxWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->func_cycle(); break;
        case 1: _t->graph_cycle(); break;
        case 2: _t->s_inc(); break;
        case 3: _t->s_decr(); break;
        case 4: _t->n_inc(); break;
        case 5: _t->n_decr(); break;
        case 6: _t->p_inc(); break;
        case 7: _t->p_decr(); break;
        case 8: _t->close(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject ApproxWindow::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_ApproxWindow.data,
      qt_meta_data_ApproxWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *ApproxWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ApproxWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_ApproxWindow.stringdata0))
        return static_cast<void*>(this);
    return QWidget::qt_metacast(_clname);
}

int ApproxWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 9)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 9;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
