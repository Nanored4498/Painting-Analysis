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

	saveAct = new QAction("&Save", this);
	saveAct->setShortcuts(QKeySequence::Save);
	saveAct->setStatusTip("Save the analysis in a SVG file");
	saveAct->setEnabled(false);
	connect(saveAct, SIGNAL(triggered()), this, SLOT(save()));

	fileMenu = menuBar()->addMenu("&File");
	fileMenu->addAction(openAct);
	fileMenu->addAction(saveAct);

	QWidget *centralW = new QWidget();
	QVBoxLayout *layout = new QVBoxLayout();
	layout->setContentsMargins(8, 2, 8, 2);
	layout->setSpacing(0);
	centralW->setLayout(layout);
	setCentralWidget(centralW);

	area = new DrawingArea();
	layout->addWidget(area);
	connect(area, SIGNAL(selected(int)), this, SLOT(shiftSelection(int)));
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

	selectionBut = new QPushButton("Compute Vanish Point");
	selectionBut->setEnabled(false);
	butLayout->addWidget(selectionBut);
	connect(selectionBut, SIGNAL(clicked()), area, SLOT(selectionAction()));

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

	rBrushSli = new LabeledSlider("Brush radius", 6, 35, RBRUSH0);
	sliLayout->addWidget(rBrushSli);
	connect(rBrushSli, SIGNAL(valueChanged(int)), area, SLOT(changeRBrush(int)));
}

void Window::open() {
	QString name = QFileDialog::getOpenFileName(this, "Open Image", QDir::currentPath(), "Images (*.png *.jpg)");
	if(!name.isEmpty()) area->loadImage(name);
}

void Window::save() {
	const QString currFile = area->getFilename();
	int li = currFile.lastIndexOf(".");
	QString name = QFileDialog::getSaveFileName(this, "Save to SVG", currFile.left(li) + ".svg", "Images (*.svg)");
	if(!name.endsWith(".svg")) name += ".svg";
	if(!name.isEmpty()) area->save_svg(name.toStdString());
}

void Window::disableButs() {
	saveAct->setEnabled(true);
	sobelBut->setEnabled(true);
	getLinesBut->setEnabled(false);
	selectionBut->setEnabled(false);
	selectionBut->setText("Compute Vanish Point");
	plotedImageBut->setCurrentIndex(0);
	QModelIndex index = plotedImageBut->model()->index(1, 0);
	plotedImageBut->model()->setData(index, 0, Qt::UserRole - 1);
}

void Window::enableGetLine() {
	getLinesBut->setEnabled(true);
	selectionBut->setEnabled(false);
	selectionBut->setText("Compute Vanish Point");
	QModelIndex index = plotedImageBut->model()->index(1, 0);
	plotedImageBut->model()->setData(index, 1|32, Qt::UserRole - 1);
}

void Window::shiftSelection(int a) {
	selectionBut->setEnabled(a != UNSELECTED);
	if(a == VANISH_POINT) selectionBut->setText("Compute Vanish Point");
	else if(a == INTERSECTIONS) selectionBut->setText("Compute Intersections");
	else if(a == UNGROUP) selectionBut->setText("Ungroup");
}

void Window::keyPressEvent(QKeyEvent *event) {
	switch(event->key()) {
	case Qt::Key_Tab:
		plotedImageBut->setCurrentIndex(1 - plotedImageBut->currentIndex());
		break;
	
	default:
		break;
	}
	QMainWindow::keyPressEvent(event);
}

Window::~Window() {
	delete fileMenu;
	delete openAct;

	delete area;

	delete getLinesBut;
	delete selectionBut;
	delete sobelBut;
	delete plotedImageBut;

	delete nLinesSli;
	delete rBrushSli;
}
