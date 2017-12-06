SRC_DIR = $$TOP_PWD/src
INCLUDEPATH += $$TOP_PWD/libs
LIBS += -larmadillo
LIBS += -lboost_system -lboost_filesystem -lX11
LIBS += -ltiff -lpng -lX11
#LIBS += -larmadillo -llapack -lblas


CONFIG *= c++11
QMAKE_RPATHDIR *= $$TOP_OUT_PWD/src
COMMON_CXXFLAGS *= -std=c++11
QMAKE_CXXFLAGS *= $$COMMON_CXXFLAGS
DEFINES += FORCE_OMP_CPU

# OPEN MP
LIBS += -fopenmp
QMAKE_CXX += -fopenmp

release {
    # Remoing other O flags
    QMAKE_CXXFLAGS_RELEASE -= -O
    QMAKE_CXXFLAGS_RELEASE -= -O1
    QMAKE_CXXFLAGS_RELEASE -= -O2

    # add the desired -O3 if not present
    QMAKE_CXXFLAGS_RELEASE *= -O3
    DEFINES *= ARMA_NO_DEBUG
}

#message(default DEFINES = $$DEFINES)
