#ifndef DRAWINGELEMENTS_H
#define DRAWINGELEMENTS_H

#include <QPoint>
#include <QLine>
#include <QtMath>

#include <lines.h>

#define PDD std::pair<double, double>

class DPoint {
public:
	DPoint(int x, int y): p0(x, y) {}

	inline void update(double s, const QPoint &sp, const QPoint &dp) {
		p = (p0-sp) * s + dp;
	}

	inline double get_dist(double px, double py) {
		double dx = p.x() - px;
		double dy = p.y() - py;
		return qSqrt(dx*dx + dy*dy);
	}

	void unselect() { selected = false; }
	void select() { selected = !selected; }
	void setGroup(int g) { group = g; }

	const QPoint& get_point() const { return p; }
	const QPoint& get_point0() const { return p0; }
	bool is_selected() const { return selected; }
	int get_group() const { return group; }

private:
	QPoint p0, p;
	bool selected = false;
	int group = 0;
};

class DLine {
public:
	DLine(const PA::Line &li, int W, int H);

	inline void update(double s, const QPoint &sp, const QPoint &dp) {
		l.setP1((l0.p1()-sp) * s + dp);
		l.setP2((l0.p2()-sp) * s + dp);
	}

	inline double get_dist(double px, double py) const {
		double d = co*px + si*py;
		double r = co*l.x1() + si*l.y1();
		return qAbs(d - r);
	}

	inline PDD getIntersection(const DLine &other) {
		double x = (other.li.b - li.b) / (li.a - other.li.a);
		double y = li.a * x + li.b;
		return {x, y};
	}

	void unselect() { selected = false; }
	void select() { selected = !selected; }
	void setGroup(int g) { group = g; }

	const QLine get_line() const { return l; }
	const QLine get_line0() const { return l0; }
	bool is_selected() const { return selected; }
	int get_group() const { return group; }
	inline void getCSR(double &c, double &s, double &r) const { c = co, s = si, r = li.rho; }
	double getTheta() const { return li.theta; }
	double getRho() const { return li.rho; }
	double getA() const { return li.a; }
	double getB() const { return li.b; }

private:
	PA::Line li;
	double co, si;
	QLine l0, l;
	bool selected = false;
	int group = -1;
};

DLine pca_pdd(const std::vector<PDD> &ps, int W, int H);

#endif // DRAWINGELEMENTS_H