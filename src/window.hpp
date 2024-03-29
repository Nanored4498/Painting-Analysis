#ifndef WINDOW_H
#define WINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QKeyEvent>

#include "drawing_area.hpp"
#include "labeled_slider.hpp"

class Window: public QMainWindow {
Q_OBJECT

public:
	Window();
	~Window();

public slots:
	void enableGetLine();
	void shiftSelection(int a);
	void disableButs();
	
protected slots:
	void open();
	void save();

protected:
	void keyPressEvent(QKeyEvent *event);

private:
	QMenu *fileMenu;
	QAction *openAct;
	QAction *saveAct;

	DrawingArea *area;

	QPushButton *getLinesBut;
	QPushButton *selectionBut;
	QPushButton *sobelBut;
	QComboBox *plotedImageBut;

	LabeledSlider *magThresholdSli, *sizeThresholdSli;
	LabeledSlider *rBrushSli;
};

#endif // WINDOW_H
