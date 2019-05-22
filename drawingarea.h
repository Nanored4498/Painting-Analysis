#ifndef DRAWINGAREA_H
#define DRAWINGAREA_H

#include <QWidget>

class DrawingArea : public QWidget {
Q_OBJECT

public:
	DrawingArea(QWidget *parent = nullptr);

	bool loadImage(const QString &fileName);

signals:

public slots:

protected:
	void paintEvent(QPaintEvent *event) override;
	void resizeEvent(QResizeEvent *event) override;

	void mousePressEvent(QMouseEvent *event) override;
	void mouseReleaseEvent(QMouseEvent *event) override;
	void wheelEvent(QWheelEvent *event) override;
	void mouseMoveEvent(QMouseEvent *eventMove) override;

private:
	void resize();

	QImage image, image0;
	double scale=1.0, sx=0.0, sy=0.0;
	void clamp_sxy();

	bool buttonPressed = false;
	QPoint pressPos;
	double press_sx=0.0, press_sy=0.0;
};

#endif // DRAWINGAREA_H
