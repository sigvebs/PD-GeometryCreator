TEMPLATE  = lib
TARGET    = meshGenerator
CONFIG   += console
CONFIG   -= app_bundle
CONFIG   -= qt

include(../default.pri)

SOURCES += meshgenrator.cpp \
	mg_functions.cpp

HEADERS +=\
	meshgenrator.h \
	mg_functions.h
