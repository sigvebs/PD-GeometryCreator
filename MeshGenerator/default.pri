INCLUDEPATH += $$TOP_PWD/libs
SRC_DIR = $$TOP_PWD/src

LIBS += -larmadillo  -lboost_system -lboost_filesystem -lX11

# OPEN MP
LIBS += -fopenmp
QMAKE_CXX += -fopenmp

QMAKE_RPATHDIR += $$TOP_OUT_PWD/src
COMMON_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
