TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp\
        lib.cpp

LIBS += -larmadillo -lblas -llapack

HEADERS += \
    integrand.h \
    methods.h \
    lib.h \
    monte_carlo.h \
    unit_tests.h

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
