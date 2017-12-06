TEMPLATE  = app
TARGET    = ../meshGenerator-3D
#TARGET    = meshGenerator-3D
CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt

DEPENDPATH += ../src
INCLUDEPATH +=  ../src

LIBS +=  -lconfig++
LIBS += ../src/libmeshGenerator.a

include(../default.pri)

SOURCES += main.cpp

include(deployment.pri)
