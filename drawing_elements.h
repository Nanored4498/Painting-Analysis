#ifndef DRAWINGELEMENTS_H
#define DRAWINGELEMENTS_H

#include <QPoint>
#include <QLine>
#include <QtMath>

#include <lines.h>

class DPoint {
public:
	DPoint(int x, int y): p0(x, y) {}

	void update(double s, const QPoint &sp, const QPoint &dp) {
		p = (p0-sp) * s + dp;
	}

	QPoint get_point() const { return p; }

private:
	QPoint p0, p;
};

class DLine {
public:
	DLine(const PA::Line &li, int W, int H): li(li) {
		int x0, y0, x1, y1;
		W --, H --;
		if(li.b > H) x0 = int((H - li.b) / li.a), y0 = H;
		else if(li.b < 0) x0 = int(- li.b / li.a), y0 = 0;
		else x0 = 0, y0 = int(li.b);
		double y = li.b + W*li.a;
		if(y > H) x1 = W + int((H - y) / li.a), y1 = H;
		else if(y < 0) x1 = W - int(y / li.a), y1 = 0;
		else x1 = W, y1 = int(y);
		l0 = QLine(x0, y0, x1, y1);
		co = qCos(li.theta);
		si = qSin(li.theta);
	}

	void update(double s, const QPoint &sp, const QPoint &dp) {
		l.setP1((l0.p1()-sp) * s + dp);
		l.setP2((l0.p2()-sp) * s + dp);
	}

	double get_dist(double px, double py) const {
		double d = co*px + si*py;
		double r = co*l.x1() + si*l.y1();
		return qAbs(d - r);
	}

	void unselect() { selected = false; }
	void select() { selected = !selected; }

	QLine get_line() const { return l; }
	bool is_selected() const { return selected; }
	void getCSR(double &c, double &s, double &r) const { c = co, s = si, r = li.rho; }

private:
	PA::Line li;
	double co, si;
	QLine l0, l;
	bool selected = false;
};

#endif // DRAWINGELEMENTS_H