#ifndef WINDOW_H
#define WINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include "drawingarea.h"

class Window: public QMainWindow {
Q_OBJECT

public:
	Window();
	~Window();

public slots:
	void open();

private:
	DrawingArea *area;

	QMenu *fileMenu;
	QAction *openAct;
};

#endif // WINDOW_H
