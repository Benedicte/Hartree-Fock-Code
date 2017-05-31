TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis.cpp \
    hartree_fock.cpp \
    coulomb_functions.cpp \
    ccd.cpp \
    test.cpp

HEADERS += \
    basis.h \
    hartree_fock.h \
    coulomb_functions.h \
    ccd.h \
    test.h

LIBS += -L/usr/local/Cellar/armadillo/7.500.0/lib -larmadillo

INCLUDEPATH += /usr/local/Cellar/armadillo/7.500.0/include
