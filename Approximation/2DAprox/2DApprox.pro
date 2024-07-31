QT += widgets
QT += gui

QMAKE_CXXFLAGS += -std=c++0x -pthread 
LIBS += -pthread

HEADERS = window.h \
          threadManager.h \
          Solve.h \
          fill_matrix.h
SOURCES = main.cpp \
          window.cpp \
          renderSurface.cpp \
          threadManager.cpp \
          Solve.cpp \
          fill_matrix.cpp