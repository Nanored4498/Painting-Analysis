#include "drawing_area.h"

#include <set>

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