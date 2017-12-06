TEMPLATE  = lib
CONFIG   += staticlib
TARGET    = meshGenerator
CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt

include(../default.pri)

SOURCES += \
	mg_functions.cpp \
    meshgenerator.cpp \
    grid.cpp

HEADERS +=\
	mg_functions.h \
    meshgenerator.h \
    grid.h
