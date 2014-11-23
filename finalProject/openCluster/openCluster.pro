LIBS += -larmadillo

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

SOURCES += \
    main.cpp \
    randomgenerator.cpp \
    distance.cpp \
    object.cpp \
    solvestep.cpp \
    system.cpp

HEADERS += \
    randomgenerator.h \
    distance.h \
    object.h \
    solvestep.h \
    system.h
