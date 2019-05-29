TEMPLATE += \
	app

QT += \
	widgets \
	gui

SOURCES += \
	main.cpp \
	labeled_slider.cpp \
	window.cpp \
	drawing_area.cpp \
	drawing_elements.cpp \
	lines.cpp \
	stb_image.cpp

HEADERS += \
	window.h \
	labeled_slider.h \
	drawing_area.h \
	drawing_elements.h \
	lines.h \
	stb_image.h \
	stb_image_write.h

QMAKE_CXXFLAGS += \
	-fopenmp

LIBS += \
	-fopenmp
