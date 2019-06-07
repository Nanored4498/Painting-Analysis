#include "drawing_area.h"

#include <QPainter>
#include <QMouseEvent>
#include <QApplication>

#include <set>
#include <algorithm>

#include <stb_image.h>

const static std::vector<QColor> colors = {Qt::blue, Qt::green, Qt::magenta, Qt::yellow, QColor(240, 0, 120),
											QColor(0, 240, 120), QColor(120, 0, 240), QColor(120, 240, 0)};

DrawingArea::DrawingArea(QWidget *parent) : QWidget(parent) {

}

bool DrawingArea::loadImage(const QString &fileName) {
	QImage newImage;
	if(!newImage.load(fileName)) return false;
	int W, H, C;
	uchar* new_im = stbi_load(fileName.toStdString().data(), &W, &C, &H, 3);
	if(!new_im) return false;
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

void DrawingArea::resizeEvent(QResizeEvent *event) {
	resize();
	QWidget::resizeEvent(event);
}

void DrawingArea::paintEvent(QPaintEvent *event) {
	if(image.isNull()) return;
	QPainter painter(this);
	painter.drawImage(im_x, im_y, image);
	double size_mul = qPow(scale, 0.3);
	QPen selCol(Qt::cyan, 1.5*size_mul);
	for(const DLine *l : lines) {
		if(l->is_selected()) painter.setPen(selCol);
		else if(l->get_group() == 0) painter.setPen(QPen(colors[0], 1.2*size_mul));
		else painter.setPen(QPen(colors[1 + (l->get_group()-1) % (colors.size()-1)], 1.2*size_mul));
		painter.drawLine(l->get_line());
	}
	QPen unselColP(Qt::red, 5*size_mul), selColP(QColor(240, 120, 0), 6*size_mul);
	for(const DPoint &p : vanishPoints) {
		painter.setPen(p.is_selected() ? selColP : unselColP);
		painter.drawPoint(p.get_point());
	}
	if(horizontalLine) {
		QPen horizonCol(QColor(0, 120, 240), 1.6*size_mul);
		painter.setPen(horizonCol);
		painter.drawLine(horizontalLine->get_line());
	}
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

void DrawingArea::eraseSobel(double px, double py) {
	double s = scale_im * scale;
	int ix = int(sx + px/s), iy = int(sy + py/s);
	int dx = 0;
	double r2 = rBrush*rBrush / (s*s);
	while(dx*dx <= r2) {
		int dy = 0;
		while(dx*dx+dy*dy <= r2) {
			if(sobelIm.rect().contains(ix+dx, iy+dy) && (sobelIm.pixel(ix+dx, iy+dy) & 0xffffff))
				pa_data->no[ix+dx + (iy+dy)*pa_data->W] = -qAbs(pa_data->no[ix+dx + (iy+dy)*pa_data->W]),
				sobelIm.setPixel(ix+dx, iy+dy, 0xffffff);
			if(sobelIm.rect().contains(ix-dx, iy+dy) && (sobelIm.pixel(ix-dx, iy+dy) & 0xffffff))
				pa_data->no[ix-dx + (iy+dy)*pa_data->W] = -qAbs(pa_data->no[ix-dx + (iy+dy)*pa_data->W]),
				sobelIm.setPixel(ix-dx, iy+dy, 0xffffff);
			if(sobelIm.rect().contains(ix+dx, iy-dy) && (sobelIm.pixel(ix+dx, iy-dy) & 0xffffff))
				pa_data->no[ix+dx + (iy-dy)*pa_data->W] = -qAbs(pa_data->no[ix+dx + (iy-dy)*pa_data->W]),
				sobelIm.setPixel(ix+dx, iy-dy, 0xffffff);
			if(sobelIm.rect().contains(ix-dx, iy-dy) && (sobelIm.pixel(ix-dx, iy-dy) & 0xffffff))
				pa_data->no[ix-dx + (iy-dy)*pa_data->W] = -qAbs(pa_data->no[ix-dx + (iy-dy)*pa_data->W]),
				sobelIm.setPixel(ix-dx, iy-dy, 0xffffff);
			dy ++;
		}
		dx ++;
	}
	resize();
}

void DrawingArea::mousePressEvent(QMouseEvent *event) {
	if(image0.isNull()) return;
	double px = event->pos().x();
	double py = event->pos().y();
	/*** Select a liene ***/
	if(event->button() == Qt::LeftButton) {
		// If shift is not pressed, unselect all elements already selected
		if(!(QApplication::keyboardModifiers() & Qt::ShiftModifier)) {
			for(DPoint &p : vanishPoints) p.unselect();
			for(DLine *l : lines) l->unselect();
		}
		// Is a point already selected
		bool psel = false;
		for(DPoint &p : vanishPoints)
			if(p.is_selected())
				psel = true;
		// Is a line already selected
		bool lsel = false;
		for(DLine *l : lines)
			if(l->is_selected())
				lsel = true;
		std::set<int> groups;
		if(psel || !lsel) {
			for(DPoint &p : vanishPoints) {
				if(p.get_dist(px, py) < 5) p.select();
				if(p.is_selected()) groups.insert(p.get_group());
			}
		}
		for(DLine *l : lines) {
			if(l->get_dist(px, py) < 5) l->select();
			if(l->is_selected()) groups.insert(l->get_group());
		}
		if(!psel && groups.count(0) > 0) {
			for(DPoint &p : vanishPoints) p.unselect();
			for(DLine *l : lines)
				if(l->get_group() != 0) l->unselect();
		} else {
			groups.erase(0);
			for(DLine *l : lines)
				if((groups.count(l->get_group()) > 0) != l->is_selected())
					l->select();
			for(DPoint &p : vanishPoints)
				if((groups.count(p.get_group()) > 0) != p.is_selected())
					p.select();
		}
		bool sel = false;
		for(DLine *l : lines)
			if(l->is_selected())
				sel = true;
		if(!sel) action = UNSELECTED;
		else if(!psel && groups.count(0) > 0) action = VANISH_POINT;
		else if(groups.size() > 1) action = INTERSECTIONS;
		else action = UNGROUP;
		emit selected(action);
		update();
	/*** Add line ***/
	} else if(plotOriginal && event->button() == Qt::RightButton) {
		unsigned int previous_size = lines.size();
		auto beg = std::remove_if(lines.begin(),
								lines.end(),
								[px, py](DLine *l) { return l->get_dist(px, py) < 4; }
					);
		for(auto l = beg; l < lines.end(); l++) (*l)->setGroup(-1);
		lines.erase(beg, lines.end());
		if(lines.size() < previous_size) {
			update();
			return;
		}
		double s = scale * scale_im;
		QPoint dp((width() - scale_im*image0.width())/2, (height() - scale_im*image0.height())/2);
		QPoint sp(sx, sy);
		for(DLine &l : candidate_lines) {
			if(l.get_group() >= 0) continue;
			l.update(s, sp, dp);
			if(l.get_dist(px, py) < 4) {
				l.setGroup(0);
				lines.push_back(&l);
				update();
				break;
			}
		}
	/*** Brush or move camera ***/
	} else if(event->button() == Qt::MiddleButton || event->button() == Qt::RightButton) {
		px -= (width() - scale_im*image0.width()) / 2;
		py -= (height() - scale_im*image0.height()) / 2;
		if(px < 0 || px > image.width() || py < 0 || py > image.height()) return;
		MidButPressed = event->button() == Qt::MiddleButton;
		rightButPressed = !MidButPressed;
		pressPos = event->pos();
		press_sx = sx;
		press_sy = sy;
		if(rightButPressed && !plotOriginal) eraseSobel(px, py);
	}
}

void DrawingArea::mouseReleaseEvent(QMouseEvent *event) {
	if(event->button() == Qt::MiddleButton) MidButPressed = false;
	else if(event->button() == Qt::RightButton) rightButPressed = false;
}

void DrawingArea::mouseMoveEvent(QMouseEvent *event) {
	if(image0.isNull()) return;
	if(MidButPressed) {
		QPoint npos = event->pos();
		double s = 1.0 / scale_im / scale;
		sx = press_sx - (npos.x() - pressPos.x()) * s;
		sy = press_sy - (npos.y() - pressPos.y()) * s;
		clamp_sxy();
		resize();
	} else if(rightButPressed && !plotOriginal) {
		double px = event->pos().x() - (width() - scale_im*image0.width()) / 2;
		double py = event->pos().y() - (height() - scale_im*image0.height()) / 2;
		eraseSobel(px, py);
		QPoint dir = event->pos() - pressPos;
		double dist = qSqrt(dir.x()*dir.x() + dir.y()*dir.y());
		for(double d = rBrush; d <= dist-rBrush; d += rBrush) {
			px -= dir.x() / dist * rBrush;
			py -= dir.y() / dist * rBrush;
			eraseSobel(px, py);
		}
		pressPos = event->pos();
	}
}

void DrawingArea::computeSobel() {
	if(!im) return;
	if(pa_data) delete pa_data;
	candidate_lines.clear();
	lines.clear();
	vanishPoints.clear();
	int W = image0.width(), H = image0.height();
	pa_data = PA::applySobel(im, W, H);
	PA::save_sobel("sobel.png", pa_data);
	sobelIm.load("sobel.png");
	emit sobelComputed();
}

void DrawingArea::findLines() {
	int W = image0.width(), H = image0.height();
	candidate_lines.clear();
	std::vector<PA::Line> ls = PA::get_lines(pa_data);
	for(const PA::Line &l : ls)
		candidate_lines.emplace_back(l, W, H);
	lines.clear();
	vanishPoints.clear();
	unsigned int i = 0;
	while(i < candidate_lines.size() && i < nbLines) {
		candidate_lines[i].setGroup(0);
		lines.push_back(&candidate_lines[i++]);
	}
	computeHorizon();
	resizeLines();
	update();
}

void DrawingArea::selectionAction() {
	if(action == VANISH_POINT) {
		double s, c, r;
		double s_cc = 0, s_cs = 0, s_ss = 0, s_cr = 0, s_sr = 0;
		std::set<int> groups;
		for(const DLine *l : lines) {
			groups.insert(l->get_group());
			if(!l->is_selected()) continue;
			l->getCSR(c, s, r);
			s_cc += c*c;
			s_cs += s*c;
			s_ss += s*s;
			s_cr += c*r;
			s_sr += s*r;
		}
		double det = s_cc * s_ss - s_cs * s_cs;
		if(qAbs(det) < 1e-6) return;
		double x = s_ss * s_cr - s_cs * s_sr;
		double y = - s_cs * s_cr + s_cc * s_sr;
		x /= det;
		y /= det;
		vanishPoints.emplace_back(x, y);
		int g = 1;
		while(groups.count(g) > 0) g++;
		for(DLine *l : lines)
			if(l->is_selected()) l->setGroup(g);
		vanishPoints[vanishPoints.size()-1].setGroup(g);
		vanishPoints[vanishPoints.size()-1].select();

	} else if(action == INTERSECTIONS) {

		std::vector<std::vector<int>> gs;
		std::vector<int> inds;
		for(int i = 0; i < int(lines.size()); i++) {
			if(lines[i]->is_selected()) {
				int g = lines[i]->get_group();
				while(g >= int(gs.size())) gs.emplace_back();
				if(gs[g].empty()) inds.push_back(g);
				gs[g].push_back(i);
			}
		}
		int I = inds.size();
		for(int i : inds) {
			double mc = 0, ms = 0;
			for(int a : gs[i]) {
				double t = 2 * lines[a]->getTheta();
				mc += qCos(t);
				ms += qSin(t);
			}
			double t = qTan(qAtan2(ms, mc) / 2);
			auto fun = [this, t](int i) {
				double c, s, r;
				lines[i]->getCSR(c, s, r);
				return r / (c + s*t);
			};
			std::sort(gs[i].begin(), gs[i].end(), [&fun](int a, int b) { return fun(a) < fun(b); });
		}
		for(int i = 0; i < I; i++) {
			for(int j = i+1; j < I; j++) {
				std::vector<std::vector<PDD>> vvp;
				for(int a : gs[inds[i]]) {
					std::vector<PDD> vp;
					for(int b : gs[inds[j]]) {
						PDD p = lines[a]->getIntersection(*lines[b]);
						vp.push_back(p);
					}
					vvp.push_back(vp);
				}
				int n = vvp.size();
				int m = vvp[0].size();
				int W = image0.width(), H = image0.height();
				for(int s = 1; s <= n+m-2; s++) {
					std::vector<PDD> ps;
					for(int a = qMax(0, s-m+1); a <= qMin(n-1, s); a++) {
						PDD p = vvp[a][s-a];
						if(p.first >= 0 && p.first < W && p.second >= 0 && p.second < H)
							ps.push_back(p);
					}
					if(ps.size() > 1){
						DLine l = pca_pdd(ps, W, H);
						l.setGroup(0);
						candidate_lines.push_back(l);
						lines.push_back(&candidate_lines[candidate_lines.size()-1]);
					}
					ps.clear();
					for(int a = qMax(0, s-m+1); a <= qMin(n-1, s); a++) {
						PDD p = vvp[n-1-a][s-a];
						if(p.first >= 0 && p.first < W && p.second >= 0 && p.second < H)
							ps.push_back(p);
					}
					if(ps.size() > 1) {
						DLine l = pca_pdd(ps, W, H);
						l.setGroup(0);
						candidate_lines.push_back(l);
						lines.push_back(&candidate_lines[candidate_lines.size()-1]);
					}
				}
			}
		}

	} else if(action == UNGROUP) {

		for(DLine *l : lines)
			if(l->is_selected()) l->setGroup(0);
		int V = vanishPoints.size();
		for(int i = 0; i < V; i++) {
			if(vanishPoints[i].is_selected()) {
				if(i != V-1) vanishPoints[i] = vanishPoints[V-1];
				vanishPoints.pop_back();
				break;
			}
		}
	}
	computeHorizon();
	resizeLines();
	update();
}

void DrawingArea::selectPlot(int original) {
	plotOriginal = (original == 0);
	resize();
}

void DrawingArea::changeNbLines(int n) {
	nbLines = n;
}

void DrawingArea::changeRBrush(int r) {
	rBrush = r;
}

void DrawingArea::computeHorizon() {
	if(vanishPoints.size() < 2) {
		if(horizontalLine) delete horizontalLine;
		horizontalLine = nullptr;
	}
	std::vector<PDD> ps;
	for(const DPoint &p : vanishPoints)
		ps.emplace_back(p.get_point0().x(), p.get_point0().y());
	DLine line = pca_pdd(ps, image0.width(), image0.height());
	if(!horizontalLine) horizontalLine = new DLine(line);
	else *horizontalLine = line;
}