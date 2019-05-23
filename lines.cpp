#include "lines.h"
#include <iostream>
#include <map>
#include <algorithm>
#include <stb_image_write.h>

#define uchar unsigned char
#ifdef MODE_PI
#define CYCLE M_PI
#else
#define CYCLE (2*M_PI)
#endif

void bilateral(uchar* in, uchar* out, int size, double sigma1, double sigma2, int W, int H) {
	sigma1 *= sigma1;
	sigma2 *= sigma2;
	#pragma omp parallel for
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < H; j++) {
			int miX = std::max(0, i - size), maX = std::min(W-1, i + size);
			int miY = std::max(0, j - size), maY = std::min(H-1, j + size);
			for(int c = 0; c < 3; c++) {
				double colij = in[(i + W*j) * 3 + c];
				double sum = 0, w = 0;
				for(int x = miX; x <= maX; x ++) {
					for(int y = miY; y <= maY; y++) {
						double col = in[(x + W*y) * 3 + c];
						double dx = x - i, dy = y - j;
						double a = std::exp( - std::abs(colij - col) / sigma1 - (dx*dx + dy*dy) / sigma2);
						sum += col * a;
						w += a;
					}
				}
				out[(i + W*j) * 3 + c] = std::max(0.0, std::min(255.0, sum / w));
			}
		}
	}
}

double sobel(uchar* im, double* &norm, double* &angle, int W, int H) {
	norm = new double[W*H];
	angle = new double[W*H];
	for(int i = 0; i < W; i++)
		norm[i] = norm[(H-1)*W+i] = 0;
	for(int i = 0; i < H; i++)
		norm[i*W] = norm[W-1+i*W] = 0;
	double max_no = 0;
	#pragma omp parallel for reduction(max: max_no)
	for(int x = 1; x < W-1; x++) {
		for(int y = 1; y < H-1; y++) {
			double no = 0;
			double co = 0, si = 0;
			int pix = 3*(x + y*W);
			for(int c = 0; c < 3; c++) {
				double ax = 3 * im[pix+3+3*W+c] + 10 * im[pix+3+c] + 3 * im[pix+3-3*W+c]
							- 3 * im[pix-3+3*W+c] - 10 * im[pix-3+c] - 3 * im[pix-3-3*W+c];
				double ay = 3 * im[pix+3+3*W+c] + 10 * im[pix+3*W+c] + 3 * im[pix-3+3*W+c]
							- 3 * im[pix+3-3*W+c] - 10 * im[pix-3*W+c] - 3 * im[pix-3-3*W+c];
				double n = ax*ax+ay*ay;
				no += n;
				double ns = std::sqrt(n);
				co += ns * ax;
				si += ns * ay;
			}
			no = std::sqrt(no);
			max_no = std::max(max_no, no);
			norm[x+y*W] = no;
			angle[x+y*W] = std::atan2(si, co);
		}
	}
	return max_no;
}

const double Mss[3][3] = {{0.0925, 0.12, 0.0925}, {0.12, 0.15, 0.12}, {0.0925, 0.12, 0.0925}};
void smooth_sep(double* norm, double* &angle, int W, int H) {
	double *a2 = new double[W*H];
	#pragma omp parallel for
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < H; j++) {
			int pix = i + j*W;
			if(i == 0 || i == W-1 || j == 0 || j == H-1) {
				a2[pix] = angle[pix];
				continue;
			}
			double co = 0, si = 0;
			for(int a = 0; a < 3; a++) {
				for(int b = 0; b < 3; b++) {
					int p2 = (i+a-1) + (j+b-1)*W;
					co += Mss[a][b] * norm[p2] * std::cos(angle[p2]);
					si += Mss[a][b] * norm[p2] * std::sin(angle[p2]);
				}
			}
			a2[pix] = std::atan2(si, co);
		}
	}
	delete[] angle;
	angle = a2;
}

void save_sobel(const char *filename, double *no, double *an, int W, int H, double m) {
	uchar *res = new uchar[W*H*3];
	int count = 0;
	for(int i = 0; i < W*H; i++) {
		double t = an[i] / CYCLE;
		#ifndef MODE_PI
		t += 0.5;
		#endif
		double n = no[i] * 255.0 / m;
		if(n < 12) n = 0;
		else count ++, n = 255.0;
		res[3*i] = (std::max(0.0, 1-t*3) + std::max(0.0, t*3-2))*n;
		res[3*i+1] = std::max(0.0, 1-std::abs(3*t-1))*n;
		res[3*i+2] = std::max(0.0, 1-std::abs(3*t-2))*n;
	}
	stbi_write_png(filename, W, H, 3, res, 0);
}

struct UF {
	std::vector<unsigned int> id, w;

	UF(unsigned int n) {
		for(unsigned int i = 0; i < n; i++) {
			id.push_back(i);
			w.push_back(1);
		}
	}

	unsigned int find(unsigned int x) {
		return id[x] == x ? x : id[x] = find(id[x]);
	}

	void merge(unsigned int a, unsigned int b){
		unsigned int fa = find(a), fb = find(b);
		if(fa == fb) return;
		if(w[fb] > w[fa]) std::swap(fb, fa);
		id[fb] = fa;
		w[fa] += w[fb];
	}

	bool same(unsigned int a, unsigned int b) {
		return find(a) == find(b);
	}
};

void clean(double *n, int W, int H, double threshold, unsigned int min_size) {
	UF u(W*H);
	int pred[4] = {-W-1, -W, -W+1, -1};
	for(int i = 1; i < W-1; i++) {
		for(int j = 1; j < H; j++) {
			int pix = i + j*W;
			if(n[pix] < threshold) continue;
			for(int k = 0; k < 4; k++) {
				int p2 = pix+pred[k];
				if(n[p2] >= threshold) u.merge(pix, p2);
			}
		}
	}
	for(int i = 0; i < W*H; i++)
		if(u.w[u.find(i)] < min_size) n[i] = 0;
}

#define TEST
PA::Line pca(std::vector<int> &ps, double *no, double *an, int W) {
	double mx = 0, my = 0, sn = 0;
	for(int pix : ps) {
		int x = pix % W, y = pix / W;
		mx += no[pix] * x;
		my += no[pix] * y;
		sn += no[pix];
	}
	mx /= sn;
	my /= sn;
	double sxx = 0, sxy = 0, syy = 0;
	#ifdef TEST
	double co = 0, si = 0;
	#endif
	for(int pix : ps) {
		double x = (pix % W) - mx, y = int(pix / W) - my;
		sxx += no[pix] * x*x;
		sxy += no[pix] * x*y;
		syy += no[pix] * y*y;
		#ifdef TEST
		double t = an[pix];
		#ifdef MODE_PI
		t *= 2;
		#endif
		co += no[pix] * std::cos(t);
		si += no[pix] * std::sin(t);
		#endif
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
	#ifdef TEST
	#ifdef MODE_PI
	double theta2 = std::atan2(si, co);
	theta2 /= 2;
	co = std::cos(theta2), si = std::sin(theta2);
	#endif
	std::swap(co, si);
	if(co < 0) co=-co;
	else si=-si;
	double n = ps.size() + 1000 * (sxx + syy) / sn / (W*W);
	double alpha = 1.0 / (1.0 + std::exp(-0.1*(n-40)));
	double a = (s*alpha + (1-alpha)*si) / (c*alpha + (1-alpha)*co);
	#else
	double a = s/c;
	#endif
	double b = my - a*mx;
	return PA::Line(a, b);
}

std::vector<PA::Line> hough_transform(double* no, double* an, int W, int H, int R, int T, double threshold,
						double* &res) {
	res = new double[R*T];
	for(int i = 0; i < R*T; i++) res[i] = 0;
	double ma = 0;
	double r_step = std::sqrt(W*W + H*H) / R;
	double t_step = 1.5 * M_PI / T;
	for(int x = 0; x < W; x++) {
		for(int y = 0; y < H; y++) {
			int pix = x + y*W;
			if(no[pix] < threshold) continue;
//			double alpha = int((std::atan2(y, x) + 0.5*M_PI) / t_step) * t_step;
//			double mit = alpha - (0.3333*T + 0.5)*t_step, mat = alpha + 0.3333*T*t_step;
			for(double tt = -0.5*M_PI + 0.5*t_step; tt < M_PI; tt += t_step) {
				double c = std::cos(tt), s = std::sin(tt);
				double r = x*c + y*s;
				if(r < 0) {
					// std::cerr << tt << " " << mit << " " << mat << "\n";
					continue;
				}
				int ri = int(r / r_step);
				int ti = int((tt + M_PI/2) / t_step);
				int p = ri + ti * R;
				res[p] += 1.0 - 0.75*std::pow(std::abs(std::sin(tt - an[pix])), 0.3);
				ma = std::max(ma, res[p]);
			}
		}
	}
	std::vector<std::pair<double, int>> lines;
	for(int i = 0; i < R*T; i++)
		if(res[i] > ma*0.41)
			lines.push_back({-res[i], i});
	std::sort(lines.begin(), lines.end());
	std::vector<std::pair<int, int>> added;
	std::vector<PA::Line> ls;
	int min_d = int(0.015*0.015*(R*R + T*T));
	for(unsigned int i = 0; i < unsigned(lines.size()); i++) {
		bool add = true;
		int x = lines[i].second % R;
		int y = lines[i].second / R;
		for(auto &other : added) {
			int dx = x - other.first;
			int dy = y - other.second;
			if(dx*dx + dy*dy < min_d) {
				add = false;
				break;
			}
		}
		if(!add) continue;
		added.push_back({x, y});
		double r = (x + 0.5) * r_step;
		double t = (y + 0.5) * t_step - M_PI/2;
		ls.emplace_back(r, t, true);
		// double max_d = 0.008 * std::sqrt(W*W + H*H);
		// double co = std::cos(t), si = std::sin(t);
		// std::vector<int> ps;
		// Line l(r, t, true);
		// for(int x = 0; x < W; x++) {
		// 	int y0 = l.a * x + l.b;
		// 	int y = y0;
		// 	while(true) {
		// 		double d = std::abs(co*x + si*y - r);
		// 		if(d > max_d || y >= H || y < 0) {
		// 			if(y >= y0) {
		// 				y = y0-1;
		// 				continue;
		// 			} else break;
		// 		}
		// 		int pix = x + y*W;
		// 		if(no[pix] >= threshold && std::abs(std::sin(t - an[pix])) < 0.2)
		// 			ps.push_back(pix);
		// 		if(y >= y0) y ++;
		// 		else y --;
		// 	}
		// }
		// Line l2 = pca(ps, no, an, W);
		// ls.push_back(l2);
	}
	return ls;
}

PA::ProblemData* applySobel(uchar* im, int W, int H, int num_smooth_pass) {
	bilateral(im, im, 3, 5, 3, W, H);
	double *no, *an;
	double m = sobel(im, no, an, W, H);
	clean(no, W, H, 0.05*m, 42);
	for(int i = 0; i < num_smooth_pass; i++)
		smooth_sep(no, an, W, H);
	#ifdef MODE_PI
	for(int i = 0; i < W*H; i++)
		if(an[i] < 0) an[i] += M_PI;
	#endif

	return new PA::ProblemData(W, H, no, an, m);
}

std::vector<PA::Line> PA::get_lines(PA::ProblemData *data) {
	double *ht;
	int R = 500, T = 360;
	std::vector<PA::Line> res = hough_transform(data->no, data->an, data->W, data->H, R, T, 0.05*data->m, ht);
	delete[] ht;

	return res;
}
