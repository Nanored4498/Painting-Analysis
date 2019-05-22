#include "drawingarea.h"
#include <QPainter>
#include <QMouseEvent>
#include <QtDebug>
#include <QtMath>
#include <QApplication>

DrawingArea::DrawingArea(QWidget *parent) : QWidget(parent) {

}

bool DrawingArea::loadImage(const QString &fileName) {
	QImage newImage;
	if(!newImage.load(fileName)) return false;
	image0 = newImage;
	resize();
	return true;
}

void DrawingArea::resize() {
	if(image0.isNull()) return;
	QRect rect(int(sx), int(sy), int(image0.width()/scale), int(image0.height()/scale));
	image = image0.copy(rect).scaled(width()-30, height()-20, Qt::KeepAspectRatio);
	update();
}

void DrawingArea::clamp_sxy() {
	sx = qMin(double(image0.width())*(1.0 - 1.0/scale), qMax(0.0, sx));
	sy = qMin(double(image0.height())*(1.0 - 1.0/scale), qMax(0.0, sy));
}

void DrawingArea::resizeEvent(QResizeEvent *event) {
	resize();
	QWidget::resizeEvent(event);
}

void DrawingArea::paintEvent(QPaintEvent *event) {
	QPainter painter(this);
	painter.drawImage((width() - image.width())/2, (height() - image.height())/2, image);
}

void DrawingArea::wheelEvent(QWheelEvent *event) {
	if(image0.isNull()) return;
	double newScale = scale * qPow(1.05, double(event->delta()) / 100.0);
	newScale = qMin(15.0, qMax(1.0, newScale));
	double px = event->pos().x() - (width() - image.width())/2;
	double py = event->pos().y() - (height() - image.height())/2;
	if(px < 0 || px > image.width() || py < 0 || py > image.height()) return;
	double s = (1 - scale/newScale) * (double(image0.width()) / double(image.width())) / scale;
	sx += px*s;
	sy += py*s;
	scale = newScale;
	clamp_sxy();
	resize();
}

void DrawingArea::mousePressEvent(QMouseEvent *event) {
	if(image0.isNull()) return;
	if(event->button() != Qt::LeftButton) return;
	double px = event->pos().x() - (width() - image.width())/2;
	double py = event->pos().y() - (height() - image.height())/2;
	if(px < 0 || px > image.width() || py < 0 || py > image.height()) return;
	buttonPressed = true;
	pressPos = event->pos();
	press_sx = sx;
	press_sy = sy;
	if(QApplication::keyboardModifiers() & Qt::ShiftModifier) {
		qInfo() << "test" << endl;
	}
}

void DrawingArea::mouseReleaseEvent(QMouseEvent *event) {
	buttonPressed = false;
}

void DrawingArea::mouseMoveEvent(QMouseEvent *event) {
	if(image0.isNull()) return;
	if(!buttonPressed) return;
	QPoint npos = event->pos();
	double s = (double(image0.width()) / double(image.width())) / scale;
	sx = press_sx - (npos.x() - pressPos.x()) * s;
	sy = press_sy - (npos.y() - pressPos.y()) * s;
	clamp_sxy();
	resize();
}
