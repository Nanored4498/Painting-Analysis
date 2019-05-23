#include "drawingarea.h"
#include <QPainter>
#include <QMouseEvent>
#include <QApplication>

#include <stb_image.h>

DrawingArea::DrawingArea(QWidget *parent) : QWidget(parent) {

}

bool DrawingArea::loadImage(const QString &fileName) {
	QImage newImage;
	if(!newImage.load(fileName)) return false;
	int W, H, C;
	uchar* new_im = stbi_load(fileName.toStdString().data(), &W, &C, &H, 3);
	if(!new_im) return false;
	image0 = newImage;
	if(im) free(im);
	im = new_im;
	if(pa_data) delete pa_data;
	lines.clear();
	vanishPoints.clear();
	emit reinitialized();
	resize();
	return true;
}

void DrawingArea::resizeLines() {
	double s = scale * scale_im;
	QPoint dp((width() - scale_im*image0.width())/2, (height() - scale_im*image0.height())/2);
	QPoint sp(sx, sy);
	for(DLine &l : lines) l.update(s, sp, dp);
	for(DPoint &p : vanishPoints) p.update(s, sp, dp);
}

void DrawingArea::resize() {
	if(image0.isNull()) return;
	int W = width(), H = height();
	scale_im = qMin(double(W)/double(image0.width()), double(H)/double(image0.height()));
	double dW = double(W) / (scale_im*image0.width()), dH = double(H) / (scale_im*image0.height());
	int rx0 = int(sx), ry0 = int(sy);
	int rw = int(image0.width()/scale), rh = int(image0.height()/scale);
	int rX0 = qMax(0, int(rx0 - 0.5*(dW-1)*rw));
	int rX1 = qMin(image0.width()-1, int(rx0 + 0.5*(dW+1)*rw));
	int rY0 = qMax(0, int(ry0 - 0.5*(dH-1)*rh));
	int rY1 = qMin(image0.height()-1, int(ry0 + 0.5*(dH+1)*rh));
	QRect rect(rX0, rY0, rX1-rX0, rY1-rY0);
	image = image0.copy(rect).scaled(W, H, Qt::KeepAspectRatio);
	im_x = qAbs(dW-1) < 1e-4 ? 0 : int(0.5 * W * (1.0 - 1.0/dW) * (1.0 - double(rx0 - rX0) / (0.5*(dW-1)*rw)));
	im_y = qAbs(dH-1) < 1e-4 ? 0 : int(0.5 * H * (1.0 - 1.0/dH) * (1.0 - double(ry0 - rY0) / (0.5*(dH-1)*rh)));
	resizeLines();
	update();
}

void DrawingArea::clamp_sxy() {
	double f = (1.0 - 1.0/scale);
	sx = qMin(double(image0.width())*f, qMax(0.0, sx));
	sy = qMin(double(image0.height())*f, qMax(0.0, sy));
}

void DrawingArea::resizeEvent(QResizeEvent *event) {
	resize();
	QWidget::resizeEvent(event);
}

void DrawingArea::paintEvent(QPaintEvent *event) {
	if(image.isNull()) return;
	QPainter painter(this);
	painter.drawImage(im_x, im_y, image);
	double size_mul = qPow(scale, 0.3);
	QPen selCol(Qt::cyan, 1.5*size_mul), notSelCol(Qt::blue, 1.2*size_mul);
	for(const DLine &l : lines) {
		painter.setPen(l.is_selected() ? selCol : notSelCol);
		painter.drawLine(l.get_line());
	}
	painter.setPen(QPen(Qt::red, 5*size_mul));
	for(const DPoint &p : vanishPoints)
		painter.drawPoint(p.get_point());
}

void DrawingArea::wheelEvent(QWheelEvent *event) {
	if(image0.isNull()) return;
	double newScale = scale * qPow(1.05, double(event->delta()) / 100.0);
	newScale = qMin(15.0, qMax(1.0, newScale));
	double px = event->pos().x() - (width() - scale_im*image0.width())/2;
	double py = event->pos().y() - (height() - scale_im*image0.height())/2;
	if(px < 0 || px >= image.width() || py < 0 || py >= image.height()) return;
	double s = (1 - scale/newScale) / scale_im / scale;
	sx += px*s;
	sy += py*s;
	scale = newScale;
	clamp_sxy();
	resize();
}

void DrawingArea::mousePressEvent(QMouseEvent *event) {
	if(image0.isNull()) return;
	double px = event->pos().x();
	double py = event->pos().y();
	if(event->button() == Qt::LeftButton) {
		if(!(QApplication::keyboardModifiers() & Qt::ShiftModifier))
			for(DLine &l : lines) l.unselect();
		for(DLine &l : lines)
			if(l.get_dist(px, py) < 6) l.select();
		bool sel = false;
		for(DLine &l : lines) {
			if(l.is_selected()) {
				sel = true;
				break;
			}
		}
		emit lineSelection(sel);
		update();
	} else if(event->button() == Qt::MiddleButton) {
		px -= (width() - image.width())/2;
		py -= (height() - image.height())/2;
		if(px < 0 || px > image.width() || py < 0 || py > image.height()) return;
		buttonPressed = true;
		pressPos = event->pos();
		press_sx = sx;
		press_sy = sy;
	}
}

void DrawingArea::mouseReleaseEvent(QMouseEvent *event) {
	if(event->button() == Qt::MiddleButton) buttonPressed = false;
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

void DrawingArea::computeSobel() {
	if(!im) return;
	if(pa_data) delete pa_data;
	lines.clear();
	vanishPoints.clear();
	int W = image0.width(), H = image0.height();
	pa_data = PA::applySobel(im, W, H, 10);
	emit sobelComputed();
}

void DrawingArea::findLines() {
	int W = image0.width(), H = image0.height();
	std::vector<PA::Line> ls = PA::get_lines(pa_data);
	lines.clear();
	vanishPoints.clear();
	for(const PA::Line &l : ls)
		lines.emplace_back(l, W, H);
	resizeLines();
	update();
}

void DrawingArea::computeVanishP() {
	double s, c, r;
	double s_cc = 0, s_cs = 0, s_ss = 0, s_cr = 0, s_sr = 0;
	for(const DLine &l : lines) {
		if(!l.is_selected()) continue;
		l.getCSR(c, s, r);
		s_cc += c*c;
		s_cs += s*c;
		s_ss += s*s;
		s_cr += c*r;
		s_sr += s*r;
	}
	double det = s_cc * s_ss - s_cs * s_cs;
	if(det == 0.0) return;
	double x = s_ss * s_cr - s_cs * s_sr;
	double y = - s_cs * s_cr + s_cc * s_sr;
	x /= det;
	y /= det;
	vanishPoints.emplace_back(x, y);
	resizeLines();
	update();
}
