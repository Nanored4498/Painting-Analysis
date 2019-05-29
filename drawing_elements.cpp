#include <drawing_elements.h>

DLine::DLine(const PA::Line &li, int W, int H): li(li) {
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

DLine pca_pdd(const std::vector<PDD> &ps, int W, int H) {
	double mx = 0, my = 0;
	for(const PDD &p : ps) {
		mx += p.first;
		my += p.second;
	}
	mx /= ps.size();
	my /= ps.size();
	double sxx = 0, sxy = 0, syy = 0;
	for(const PDD &p : ps) {
		double x = p.first - mx;
		double y = p.second - my;
		sxx += x*x;
		sxy += x*y;
		syy += y*y;
	}
	double ta2 = 2 * sxy / (sxx - syy);
	double theta = std::atan(ta2) / 2;
	double c = std::cos(theta), s = std::sin(theta);
	double l1 = (sxx * c + sxy * s) / c;
	if(2*l1 < sxx + syy) {
		std::swap(s, c);
		if(c < 0) c=-c;
		else s=-s;
	}
	double a = s/c;
	double b = my - a*mx;
	return DLine(PA::Line(a, b), W, H);
}