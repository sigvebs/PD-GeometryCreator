include(../default.pri)
TEMPLATE  = app
TARGET    = meshGenerator
CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt

DEPENDPATH += . ../src
INCLUDEPATH +=  ../src

INCLUDEPATH += $$SRC_DIR
SOURCES += main.cpp
LIBS+=  -lconfig++ -L../src/ -lmeshGenerator
#LIBS += -L$$SRC_DIR/src -lmeshGenerator
message(mdef = $$SRC_DIR)
