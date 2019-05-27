#include <window.h>

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
	layout->setContentsMargins(8, 2, 8, 2);
	layout->setSpacing(0);
	centralW->setLayout(layout);
	setCentralWidget(centralW);

	area = new DrawingArea();
	layout->addWidget(area);
	connect(area, SIGNAL(lineSelection(bool)), this, SLOT(enableVanishPoint(bool)));
	connect(area, SIGNAL(reinitialized()), this, SLOT(disableButs()));

	QWidget *butsW = new QWidget();
	butsW->setFixedHeight(42);
	QHBoxLayout *butLayout = new QHBoxLayout();
	butLayout->setContentsMargins(0, 0, 0, 0);
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

	plotedImageBut = new QComboBox();
	plotedImageBut->addItem("Original");
	plotedImageBut->addItem("Sobel");
	QModelIndex index = plotedImageBut->model()->index(1, 0);
	plotedImageBut->model()->setData(index, 0, Qt::UserRole - 1);
	butLayout->addWidget(plotedImageBut);
	connect(plotedImageBut, SIGNAL(currentIndexChanged(int)), area, SLOT(selectPlot(int)));

	QWidget *slidW = new QWidget();
	slidW->setFixedHeight(40);
	QHBoxLayout *sliLayout = new QHBoxLayout();
	sliLayout->setContentsMargins(0, 0, 0, 0);
	slidW->setLayout(sliLayout);
	layout->addWidget(slidW);

	nLinesSli = new LabeledSlider("Number of lines", 5, 40, NBLINES0);
	sliLayout->addWidget(nLinesSli);
	connect(nLinesSli, SIGNAL(valueChanged(int)), area, SLOT(changeNbLines(int)));

	rBrushSli = new LabeledSlider("Brush radius", 2, 30, RBRUSH0);
	sliLayout->addWidget(rBrushSli);
	connect(rBrushSli, SIGNAL(valueChanged(int)), area, SLOT(changeRBrush(int)));
}

void Window::open() {
	QString name = QFileDialog::getOpenFileName(this, "Open Image", QDir::currentPath(), "Images (*.png *.jpg)");
	if(!name.isEmpty()) area->loadImage(name);
}

void Window::disableButs() {
	sobelBut->setEnabled(true);
	getLinesBut->setEnabled(false);
	vanishPointBut->setEnabled(false);
	plotedImageBut->setCurrentIndex(0);
	QModelIndex index = plotedImageBut->model()->index(1, 0);
	plotedImageBut->model()->setData(index, 0, Qt::UserRole - 1);
}

void Window::enableGetLine() {
	getLinesBut->setEnabled(true);
	vanishPointBut->setEnabled(false);
	QModelIndex index = plotedImageBut->model()->index(1, 0);
	plotedImageBut->model()->setData(index, 1|32, Qt::UserRole - 1);
}

void Window::enableVanishPoint(bool enabled) {
	vanishPointBut->setEnabled(enabled);
}

Window::~Window() {
	delete fileMenu;
	delete openAct;

	delete area;

	delete getLinesBut;
	delete vanishPointBut;
	delete sobelBut;
	delete plotedImageBut;

	delete nLinesSli;
	delete rBrushSli;
}
