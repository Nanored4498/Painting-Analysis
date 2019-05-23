TEMPLATE += \
    app

QT += \
    widgets \
    gui

SOURCES += \
    main.cpp \
    window.cpp \
    drawingarea.cpp \
    lines.cpp \
    stb_image.cpp

HEADERS += \
    window.h \
    drawingarea.h \
    lines.h \
    stb_image.h \
    stb_image_write.h

QMAKE_CXXFLAGS += \
    -fopenmp

LIBS += \
    -fopenmp
