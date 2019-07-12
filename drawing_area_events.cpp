#include "drawing_area.h"

#include <QPainter>
#include <QMouseEvent>
#include <QApplication>

#include <set>

void DrawingArea::resizeEvent(QResizeEvent *event) {
	resize();
	QWidget::resizeEvent(event);
}

void DrawingArea::paintEvent(__attribute__((unused)) QPaintEvent *event) {
	if(image.isNull()) return;
	QPainter painter(this);
	painter.drawImage(im_x, im_y, image);
	double size_mul = qPow(scale, 0.3);
	/*** Sobel zone ***/
	if(pa_data == nullptr) {
		painter.setPen(QPen(Qt::white, 1.2*size_mul));
		for(int i = 1; i < int(zonePoints.size()); i++)
			painter.drawLine(zonePoints[i-1].get_point(), zonePoints[i].get_point());
		if(zonePoints.size() >= 2)
			painter.drawLine(zonePoints[0].get_point(), zonePoints.back().get_point());
	}
	/*** Parallel lines ***/
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
	double px = event->pos().x();
	double py = event->pos().y();
	if(px < im_x || px >= im_x+image.width() || py < im_y || py >= im_y+image.height()) return;
	px -= (width() - scale_im*image0.width()) / 2.0;
	py -= (height() - scale_im*image0.height()) / 2.0;
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
	/*** Select a line ***/
	if(event->button() == Qt::LeftButton) {
		if(QApplication::keyboardModifiers() & Qt::ControlModifier) {
			if(!pa_ldata) return;
			double dx = (width() - scale_im*image0.width())/2.0 , dy = (height() - scale_im*image0.height())/2.0;
			double s = scale * scale_im;
			double x0 = sx + (px - dx) / s, y0 = sy + (py - dy) / s;
			std::vector<std::pair<double, std::pair<double, double>>> goods;
			for(int t = 0; t < pa_ldata->T; t++) {
				double tt = 1.5 * M_PI * (t+0.5) / pa_ldata->T - 0.5 * M_PI;
				double c = qCos(tt), s = qSin(tt);
				double rr = x0 * c + y0 * s;
				int r = rr / qSqrt(image0.width()*image0.width() + image0.height()*image0.height()) * pa_ldata->R;
				if(r < 0 || r >= pa_ldata->R) continue;
				qInfo("%f %f %f", tt, rr, pa_ldata->score[r + pa_ldata->R*t]);
				goods.push_back({pa_ldata->score[r + pa_ldata->R*t], {tt, rr}});
			}
			std::sort(goods.begin(), goods.end());
			for(int i = 0; i < qMin(20, (int) goods.size()); i++) {
				double tt = goods[i].second.first, rr = goods[i].second.second;
				DLine *l = new DLine(PA::Line(rr, tt, true), pa_data->W, pa_data->H);
				l->setGroup(42);
				lines.push_back(l);
			}
			update();
			return;
		}
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
	/*** Add line or select Sobel zone ***/
	} else if(plotOriginal && event->button() == Qt::RightButton) {
		/*** Select Soble Zone ***/
		if(pa_data == nullptr) {
			if(px < 0.015 * width()) px = 0;
			if(py < 0.015 * height()) py = 0;
			if(px > 0.985 * width()) px = width();
			if(py > 0.985 * height()) py = height();
			double dx = (width() - scale_im*image0.width())/2.0 , dy = (height() - scale_im*image0.height())/2.0;
			double s = scale * scale_im;
			zonePoints.emplace_back(sx + (px - dx) / s, sy + (py - dy) / s);
			resizeLines();
			update();
			return;
		}
		/*** Add line ****/
		press_sx = px, press_sy = py;
	/*** Brush or move camera ***/
	} else if(event->button() == Qt::MiddleButton || event->button() == Qt::RightButton) {
		if(px < im_x || px >= im_x+image.width() || py < im_y || py >= im_y+image.height()) return;
		px -= (width() - scale_im*image0.width()) / 2;
		py -= (height() - scale_im*image0.height()) / 2;
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
	else if(event->button() == Qt::RightButton) {
		rightButPressed = false;
		if(plotOriginal && pa_data != nullptr) {
			double px = event->pos().x();
			double py = event->pos().y();
			double s = scale * scale_im;
			unsigned int previous_size = lines.size();
			for(int i = 0; i < (int) lines.size(); i++) {
				if(lines[i]->get_dist(px, py) + lines[i]->get_dist(press_sx, press_sy) < 8) {
					int g = lines[i]->get_group();
					lines[i]->setGroup(-1);
					lines[i] = lines.back();
					lines.pop_back();
					// If it remains only one line in the group then we delete the group
					if(g > 0) {
						int j = -1;
						for(int k = 0; k < (int) lines.size(); k++) {
							if(lines[k]->get_group() != g) continue;
							if(j < 0) j = k;
							else { j = -1; break; }
						}
						if(j >= 0) lines[j]->setGroup(0);
					}
					/******************************************************************/
					i--;
				}
			}
			if(lines.size() < previous_size) {
				update();
				return;
			}
			QPoint dp((width() - scale_im*image0.width())/2.0, (height() - scale_im*image0.height())/2.0);
			QPoint sp(sx, sy);
			for(DLine &l : candidate_lines) {
				if(l.get_group() >= 0) continue;
				l.update(s, sp, dp);
				if(l.get_dist(px, py) + l.get_dist(press_sx, press_sy) < 8) {
					l.setGroup(0);
					lines.push_back(&l);
					update();
					break;
				}
			}
		}
	}
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
