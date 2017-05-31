TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    basis.cpp \
    coulomb_functions.cpp \
    hartree_fock_equations.cpp

HEADERS += \
    basis.h \
    coulomb_functions.h \
    hartree_fock_equations.h


LIBS += -L/usr/local/Cellar/armadillo/7.500.0/lib -larmadillo

INCLUDEPATH += /usr/local/Cellar/armadillo/7.500.0/include
