#include "window.h"
#include <QApplication>
#include <QMenuBar>
#include <QFileDialog>

Window::Window(): QMainWindow() {
	setWindowTitle("Painting Analysis");
	setMinimumSize(100, 100);
	resize(400, 300);

	area = new DrawingArea();
	setCentralWidget(area);

	openAct = new QAction("&Open", this);
	openAct->setShortcuts(QKeySequence::Open);
	openAct->setStatusTip("Open an image");
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

	fileMenu = menuBar()->addMenu("&File");
	fileMenu->addAction(openAct);
}

void Window::open() {
	QString name = QFileDialog::getOpenFileName(this, "Open Image", QDir::currentPath(), "Images (*.png *.jpg)");
	if(!name.isEmpty()) area->loadImage(name);
}

Window::~Window() {
	delete area;
	delete fileMenu;
	delete openAct;
}
