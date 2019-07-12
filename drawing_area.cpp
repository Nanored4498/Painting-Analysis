#include "drawing_area.h"

#include <QApplication>

#include "stb_image.h"
#include <fstream>

const std::vector<QColor> DrawingArea::colors = {
	Qt::blue, Qt::green, Qt::magenta, Qt::yellow, QColor(240, 0, 120),
	QColor(0, 240, 120), QColor(120, 0, 240), QColor(120, 240, 0)
};

DrawingArea::DrawingArea(QWidget *parent) : QWidget(parent) {

}

bool DrawingArea::loadImage(const QString &fileName) {
	QImage newImage;
	if(!newImage.load(fileName)) {
		qDebug("\033[31mError while loading the file: %s\033[0m", fileName.toStdString().c_str());
		return false;
	}
	int W, H, C;
	uchar* new_im = stbi_load(fileName.toStdString().data(), &W, &C, &H, 3);
	if(!new_im) return false;
	filename0 = fileName;
	image0 = newImage;
	if(im) {
		free(im);
		im = nullptr;
	}
	im = new_im;
	if(pa_data) {
		delete pa_data;
		pa_data = nullptr;
	}
	zonePoints.clear();
	candidate_lines.clear();
	lines.clear();
	vanishPoints.clear();
	computeHorizon();
	emit reinitialized();
	scale = 1.0;
	sx = 0.0;
	sy = 0.0;
	plotOriginal = true;
	resize();
	return true;
}

void DrawingArea::resizeLines() {
	double s = scale * scale_im;
	QPoint dp((width() - scale_im*image0.width())/2, (height() - scale_im*image0.height())/2);
	QPoint sp(sx, sy);
	for(DPoint &p : zonePoints) p.update(s, sp, dp);
	for(DLine *l : lines) l->update(s, sp, dp);
	for(DPoint &p : vanishPoints) p.update(s, sp, dp);
	if(horizontalLine) horizontalLine->update(s, sp, dp);
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
	if(plotOriginal) image = image0.copy(rect).scaled(W, H, Qt::KeepAspectRatio);
	else image = sobelIm.copy(rect).scaled(W, H, Qt::KeepAspectRatio);
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

bool DrawingArea::isSobelPixel(int x, int y) {
	return sobelIm.rect().contains(x, y) && qAbs(pa_data->no[x + y*pa_data->W]) > 1e-5;
}

void DrawingArea::shiftSobelPixel(int x, int y, bool unbrush) {
	if(unbrush) {
		pa_data->no[x + y*pa_data->W] = qAbs(pa_data->no[x + y*pa_data->W]),
		sobelIm.setPixel(x, y, PA::pixelColor(x, y, pa_data));
	} else {
		pa_data->no[x + y*pa_data->W] = -qAbs(pa_data->no[x + y*pa_data->W]),
		sobelIm.setPixel(x, y, 0xffffff);
	}
}

void DrawingArea::eraseSobel(double px, double py) {
	double s = scale_im * scale;
	int ix = int(sx + px/s), iy = int(sy + py/s);
	int dx = 0;
	double r2 = rBrush*rBrush / (s*s);
	bool unbrush = QApplication::keyboardModifiers() & Qt::ShiftModifier;
	while(dx*dx <= r2) {
		int dy = 0;
		while(dx*dx+dy*dy <= r2) {
			if(isSobelPixel(ix+dx, iy+dy)) shiftSobelPixel(ix+dx, iy+dy, unbrush);
			if(isSobelPixel(ix-dx, iy+dy)) shiftSobelPixel(ix-dx, iy+dy, unbrush);
			if(isSobelPixel(ix+dx, iy-dy)) shiftSobelPixel(ix+dx, iy-dy, unbrush);
			if(isSobelPixel(ix-dx, iy-dy)) shiftSobelPixel(ix-dx, iy-dy, unbrush);
			dy ++;
		}
		dx ++;
	}
	resize();
}

void DrawingArea::computeHorizon() {
	if(vanishPoints.size() < 2) {
		if(horizontalLine) delete horizontalLine;
		horizontalLine = nullptr;
		return;
	}
	std::vector<PDD> ps;
	for(const DPoint &p : vanishPoints)
		ps.emplace_back(p.get_point0().x(), p.get_point0().y());
	DLine line = pca_pdd(ps, image0.width(), image0.height());
	if(!horizontalLine) horizontalLine = new DLine(line);
	else *horizontalLine = line;
}

void DrawingArea::save_svg(const std::string &filename) const {
	double xmin = 0, xmax = image0.width(), ymin = 0, ymax = image0.height();
	for(const DPoint &p : vanishPoints) {
		xmin = qMin((double) p.get_point0().x(), xmin);
		xmax = qMax((double) p.get_point0().x(), xmax);
		ymin = qMin((double) p.get_point0().y(), ymin);
		ymax = qMax((double) p.get_point0().y(), ymax);
	}
	double dw = (xmax - xmin) * 0.05, dh = (ymax - ymin) * 0.05;
	xmin -= dw, xmax += dw, ymin -= dh, ymax += dh;
	std::ofstream res;
	res.open(filename);
	res << "<svg width=\"" << xmax-xmin << "\" height=\"" << ymax-ymin << "\" version=\"1.1\""
		<< " xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";
	xmin -= 0.5, ymin -= 0.5;
	double diag = qSqrt(qPow(xmax-xmin, 2.0) + qPow(ymax-ymin, 2.0));
	double line_width = 1.14 + diag * 19e-5;
	double horizon_width = 1.88 + diag * 30e-5;
	double point_width = 2.67 + diag * 44e-5;
	res << "\t<image xlink:href=\"" << filename0.toStdString() << "\" x=\"" << int(-xmin) << "\" y=\"" << int(-ymin) << "\""
		<< " height=\"" << image0.height() << "\" width=\"" << image0.width() << "\" />\n";
	for(const DLine* l : lines) {
		const QLine l0 = l->get_line0();
		QColor col = l->get_group() == 0 ? colors[0] : colors[1 + (l->get_group() - 1) % (colors.size()-1)];
		res << "\t<line x1=\"" << int(l0.x1()-xmin) << "\" y1=\"" << int(l0.y1()-ymin) << "\""
			<< " x2=\"" << int(l0.x2()-xmin) << "\" y2=\"" << int(l0.y2()-ymin)
			<< "\" style=\"stroke:rgb(" << col.red() << ", " << col.green() << ", " << col.blue() << "); stroke-width:" << line_width << "\" />\n";
	}
	if(horizontalLine != nullptr) {
		const QLine l0 = horizontalLine->get_line0();
		res << "\t<line x1=\"" << int(l0.x1()-xmin) << "\" y1=\"" << int(l0.y1()-ymin) << "\""
			<< " x2=\"" << int(l0.x2()-xmin) << "\" y2=\"" << int(l0.y2()-ymin)
			<< "\" style=\"stroke:rgb(0, 120, 240); stroke-width:" << horizon_width << "\" />\n";
	}
	for(const DPoint &p : vanishPoints) {
		const QPoint qp = p.get_point0();
		res << "\t<circle cx=\"" << int(qp.x()-xmin) << "\" cy=\"" << int(qp.y()-ymin) << "\" r=\"" << point_width << "\" fill=\"red\" />\n";
	}
	res << "</svg>\n";
	res.close();
} 