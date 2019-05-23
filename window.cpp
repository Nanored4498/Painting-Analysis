#include "window.h"
#include <QApplication>
#include <QMenuBar>
#include <QFileDialog>
#include <QVBoxLayout>
#include <QHBoxLayout>

Window::Window(): QMainWindow() {
	setWindowTitle("Painting Analysis");
	setMinimumSize(300, 150);
	resize(400, 300);

	openAct = new QAction("&Open", this);
	openAct->setShortcuts(QKeySequence::Open);
	openAct->setStatusTip("Open an image");
	connect(openAct, SIGNAL(triggered()), this, SLOT(open()));

	fileMenu = menuBar()->addMenu("&File");
	fileMenu->addAction(openAct);

	QWidget *centralW = new QWidget();
	QVBoxLayout *layout = new QVBoxLayout();
	centralW->setLayout(layout);
	setCentralWidget(centralW);

	area = new DrawingArea();
	layout->addWidget(area);
	connect(area, SIGNAL(lineSelection(bool)), this, SLOT(enableVanishPoint(bool)));
	connect(area, SIGNAL(reinitialized()), this, SLOT(disableButs()));

	QWidget *butsW = new QWidget();
	butsW->setFixedHeight(42);
	QHBoxLayout *butLayout = new QHBoxLayout();
	butsW->setLayout(butLayout);
	layout->addWidget(butsW);

	sobelBut = new QPushButton("Compute Sobel");
	sobelBut->setEnabled(false);
	butLayout->addWidget(sobelBut);
	connect(sobelBut, SIGNAL(clicked()), area, SLOT(computeSobel()));
	connect(area, SIGNAL(sobelComputed()), this, SLOT(enableGetLine()));

	getLinesBut = new QPushButton("Find Lines");
	getLinesBut->setEnabled(false);
	butLayout->addWidget(getLinesBut);
	connect(getLinesBut, SIGNAL(clicked()), area, SLOT(findLines()));

	vanishPointBut = new QPushButton("Compute Vanish Point");
	vanishPointBut->setEnabled(false);
	butLayout->addWidget(vanishPointBut);
	connect(vanishPointBut, SIGNAL(clicked()), area, SLOT(computeVanishP()));
}

void Window::open() {
	QString name = QFileDialog::getOpenFileName(this, "Open Image", QDir::currentPath(), "Images (*.png *.jpg)");
	if(!name.isEmpty()) area->loadImage(name);
}

void Window::disableButs() {
	sobelBut->setEnabled(true);
	getLinesBut->setEnabled(false);
	vanishPointBut->setEnabled(false);
}

void Window::enableGetLine() {
	getLinesBut->setEnabled(true);
	vanishPointBut->setEnabled(false);
}

void Window::enableVanishPoint(bool enabled) {
	vanishPointBut->setEnabled(enabled);
}

Window::~Window() {
	delete fileMenu;
	delete openAct;
	delete area;
	delete getLinesBut;
}
