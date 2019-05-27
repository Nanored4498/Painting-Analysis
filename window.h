#ifndef WINDOW_H
#define WINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>

#include <drawing_area.h>
#include <labeled_slider.h>

class Window: public QMainWindow {
Q_OBJECT

public:
	Window();
	~Window();

public slots:
	void open();
	void enableGetLine();
	void enableVanishPoint(bool enabled);
	void disableButs();

private:
	QMenu *fileMenu;
	QAction *openAct;

	DrawingArea *area;

	QPushButton *getLinesBut;
	QPushButton *vanishPointBut;
	QPushButton *sobelBut;
	QComboBox *plotedImageBut;

	LabeledSlider *nLinesSli;
	LabeledSlider *rBrushSli;
};

#endif // WINDOW_H
