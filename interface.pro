TEMPLATE += \
	app

QT += \
	widgets \
	gui

SOURCES += \
	src/main.cpp \
	src/labeled_slider.cpp \
	src/window.cpp \
	src/drawing_area.cpp \
	src/drawing_area_events.cpp \
	src/drawing_area_slots.cpp \
	src/drawing_elements.cpp \
	src/lines.cpp \
	src/stb_image/stb_image.cpp

HEADERS += \
	src/window.hpp \
	src/labeled_slider.hpp \
	src/drawing_area.hpp \
	src/drawing_elements.hpp \
	src/lines.hpp \
	src/stb_image/stb_image.h \
	src/stb_image/stb_image_write.h

QMAKE_CXXFLAGS += \
	-fopenmp

LIBS += \
	-fopenmp

CONFIG(release, debug|release) {
    CONFIG += optimize_full
}