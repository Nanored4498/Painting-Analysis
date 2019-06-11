#ifndef DRAWINGAREA_H
#define DRAWINGAREA_H

#include <QWidget>
#include <QtMath>
#include <QtDebug>

#include <vector>

#include <drawing_elements.h>

#define NBLINES0 20
#define RBRUSH0 8

#define UNSELECTED 0
#define VANISH_POINT 1
#define INTERSECTIONS 2
#define UNGROUP 3

class DrawingArea : public QWidget {
Q_OBJECT

public:
	DrawingArea(QWidget *parent = nullptr);

	bool loadImage(const QString &fileName);

signals:
	void reinitialized();
	void sobelComputed();
	void selected(int a);

public slots:
	void findLines();
	void computeSobel();
	void selectionAction();
	void selectPlot(int original);
	void changeNbLines(int n);
	void changeRBrush(int r);

protected:
	void paintEvent(QPaintEvent *event) override;
	void resizeEvent(QResizeEvent *event) override;

	void mousePressEvent(QMouseEvent *event) override;
	void mouseReleaseEvent(QMouseEvent *event) override;
	void wheelEvent(QWheelEvent *event) override;
	void mouseMoveEvent(QMouseEvent *eventMove) override;

private:
	void resize();
	void resizeLines();
	void shiftSobelPixel(int x, int y, bool unbrush=false);
	void eraseSobel(double px, double py);
	void computeHorizon();

	QImage image, image0;
	uchar* im = nullptr;
	double scale, sx, sy;
	int im_x, im_y;
	double scale_im;
	void clamp_sxy();

	PA::ProblemData *pa_data = nullptr;
	QImage sobelIm;
	bool plotOriginal = true;

	bool MidButPressed = false;
	QPoint pressPos;
	double press_sx=0.0, press_sy=0.0;

	double rBrush = RBRUSH0;
	bool rightButPressed = false;

	int action = 0;

	unsigned int nbLines = NBLINES0;
	std::vector<DLine> candidate_lines;
	std::vector<DLine*> lines;
	std::vector<DPoint> vanishPoints;
	DLine* horizontalLine = nullptr;
};

#endif // DRAWINGAREA_H
