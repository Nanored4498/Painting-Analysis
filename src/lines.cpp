#include "lines.hpp"
#include <iostream>
#include <map>
#include <algorithm>
#include "stb_image/stb_image_write.h"
#include <chrono>

#define uchar unsigned char
#ifdef MODE_PI
#define CYCLE M_PI
#else
#define CYCLE (2*M_PI)
#endif

// #define DEBUG_BIL
// #define DEBUG_SOB
#define DEBUG_CLE
// #define DEBUG_HOU

void bilateral(Color* in, Color* out, int size, double sigma_col, double sigma_dis, int W, int H, bool *mask=nullptr) {
	sigma_col *= sigma_col;
	sigma_dis *= sigma_dis;
	#pragma omp parallel for
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < H; j++) {
			int pix = i + W*j;
			Color new_col(0, 0, 0);
			if(mask != nullptr && !mask[pix]) {
				out[pix] = new_col;
				continue;
			}
			int miX = std::max(0, i - size), maX = std::min(W-1, i + size);
			int miY = std::max(0, j - size), maY = std::min(H-1, j + size);
			double w = 0;
			for(int x = miX; x <= maX; x ++) {
				for(int y = miY; y <= maY; y++) {
					double dc2 = 0;
					int pix2 = x + W*y;
					for(int c = 0; c < 3; c++) {
						double dc = in[pix2].get(c) - in[pix].get(c);
						dc2 += dc * dc;
					}
					double dx = x - i, dy = y - j;
					double a = std::exp( - dc2 / sigma_col - (dx*dx + dy*dy) / sigma_dis);
					w += a;
					new_col += in[pix2] * a;
				}
			}
			out[pix] = new_col / w;
		}
	}
}

double sobel(Color* im, double* &norm, double* &angle, int W, int H, bool *mask=nullptr) {
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
			int pix = x + y*W;
			if(mask != nullptr) {
				bool bad = false;
				for(int i = -1; i <= 1; i++)
					for(int j = -W; j <= W; j += W)
						if(!mask[pix+i+j]) bad = true;
				if(bad) {
					norm[pix] = angle[pix] = 0;
					continue;
				}
			}
			double no = 0;
			double co = 0, si = 0;
			Color cax = 3 * im[pix+1+W] + 10 * im[pix+1] + 3 * im[pix+1-W]
						- 3 * im[pix-1+W] - 10 * im[pix-1] - 3 * im[pix-1-W];
			Color cay = 3 * im[pix+1+W] + 10 * im[pix+W] + 3 * im[pix-1+W]
						- 3 * im[pix+1-W] - 10 * im[pix-W] - 3 * im[pix-1-W];
			for(int c = 0; c < 3; c++) {
				double ax = cax.get(c), ay = cay.get(c);
				#ifdef MODE_PI
				if(ax < 0) ax = -ax, ay = -ay;
				#endif
				double n = ax*ax+ay*ay;
				no += n;
				double ns = std::sqrt(n);
				co += ns * ax;
				si += ns * ay;
			}
			no = std::sqrt(no);
			max_no = std::max(max_no, no);
			norm[pix] = no;
			angle[pix] = std::atan2(si, co);
		}
	}
	return max_no;
}

const double Mss[3][3] = {{0.0925, 0.12, 0.0925}, {0.12, 0.15, 0.12}, {0.0925, 0.12, 0.0925}};
void smooth_sep(double* norm, double* angle, int W, int H, double* res, double m) {
	#pragma omp parallel for
	for(int i = 0; i < W; i++) {
		for(int j = 0; j < H; j++) {
			int pix = i + j*W;
			if(norm[pix] < 0.001*m) continue;
			if(i == 0 || i == W-1 || j == 0 || j == H-1) {
				res[pix] = angle[pix];
				continue;
			}
			double co = 0, si = 0;
			for(int a = 0; a < 3; a++) {
				for(int b = 0; b < 3; b++) {
					int p2 = (i+a-1) + (j+b-1)*W;
					if(norm[p2] < 0.001*m) continue;
					#ifdef MODE_PI
					co += Mss[a][b] * norm[p2] * std::cos(2.0*angle[p2]);
					si += Mss[a][b] * norm[p2] * std::sin(2.0*angle[p2]);
					#else
					co += Mss[a][b] * norm[p2] * std::cos(angle[p2]);
					si += Mss[a][b] * norm[p2] * std::sin(angle[p2]);
					#endif
				}
			}
			res[pix] = std::atan2(si, co);
			#ifdef MODE_PI
			res[pix] /= 2.0;
			#endif
		}
	}
}


uint pixelColori(int i, PA::ProblemData *data) {
	double t = data->an[i] / CYCLE;
	#ifndef MODE_PI
	t += 0.5;
	#endif
	double n = data->no[i] / data->m;
	if(n > 0) n = (255 - MIN_SOBEL_INTENSITY) * std::pow(n, 0.5) + MIN_SOBEL_INTENSITY;
	uint res = (std::max(0.0, 1-t*3) + std::max(0.0, t*3-2))*n;
	res = (res << 8) + std::max(0.0, 1-std::abs(3*t-1))*n;
	res = (res << 8) + std::max(0.0, 1-std::abs(3*t-2))*n;
	return res;
}

uint PA::pixelColor(int x, int y, PA::ProblemData *data) {
	return pixelColori(x + y * data->W, data);
}

void PA::save_sobel(const char *filename, PA::ProblemData *data) {
	uchar *res = new uchar[data->W*data->H*3];
	for(int i = 0; i < data->W*data->H; i++) {
		uint col = pixelColori(i, data);
		res[3*i+2] = col & 0xff;
		res[3*i+1] = (col >>= 8) & 0xff;
		res[3*i] = (col >>= 8) & 0xff;
	}
	stbi_write_png(filename, data->W, data->H, 3, res, 0);
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

PA::ProblemData* PA::applySobelToBil(Color* im, int W, int H, bool *mask, double threshold, double size_threshold) {
	auto time = std::chrono::high_resolution_clock::now();
	double diag = std::sqrt(W*W + H*H);
	double sigma_col = 250.0 / pow(diag, 0.4);
	int num_smooth_pass = 0.05 * std::pow(diag, 0.63);
	double *no, *an;
	
	double m = sobel(im, no, an, W, H, mask);
	PA::ProblemData *res = new PA::ProblemData(W, H, no, an, im, m);
	auto time2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dtime = time2 - time;
	std::cerr << "Sobel filter: " << dtime.count() << std::endl;
	time = std::chrono::high_resolution_clock::now();
	/* an in (-pi/2, pi/2) */

	#if defined DEBUG_SOB || defined DEBUG_CLE
	#ifdef MODE_PI
	for(int i = 0; i < W*H; i++)
		if(an[i] < 0) an[i] += M_PI;
	#endif
	#endif
	#ifdef DEBUG_SOB
	save_sobel("sobel0.png", res);
	#endif

	clean(no, W, H, 0.001*threshold*sigma_col*m, 0.01*size_threshold*diag);
	time2 = std::chrono::high_resolution_clock::now();
	dtime = time2 - time;
	std::cerr << "Cleaning: " << dtime.count() << std::endl;
	time = std::chrono::high_resolution_clock::now();

	#ifdef DEBUG_CLE
	save_sobel("sobel1.png", res);
	#endif

	double *an2 = new double[W*H];
	for(int i = 0; i < num_smooth_pass; i++) {
		smooth_sep(no, an, W, H, an2, m);
		std::swap(an, an2);
	}
	res->an = an;
	delete[] an2;
	time2 = std::chrono::high_resolution_clock::now();
	dtime = time2 - time;
	std::cerr << "Smoothing (" << num_smooth_pass << "): " << dtime.count() << std::endl;
	time = std::chrono::high_resolution_clock::now();

	#ifdef MODE_PI
	for(int i = 0; i < W*H; i++)
		if(an[i] < 0) an[i] += M_PI;
	#endif

	return res;
}

PA::ProblemData* PA::applySobel(uchar* im, int W, int H, bool *mask, double threshold, double size_threshold) {
	double diag = std::sqrt(W*W + H*H);
	double sigma = std::pow(diag, 0.3) * 0.3;
	double sigma_col = 250.0 / pow(diag, 0.4);
	auto time = std::chrono::high_resolution_clock::now();
	
	Color *im2 = new Color[W*H];
	#pragma omp parallel for
	for(int i = 0; i < W*H; i++)
		if(mask == nullptr || mask[i]) im2[i] = rgb_to_lab(Color(im[3*i], im[3*i+1], im[3*i+2]));
		else im2[i] = Color(0, 0, 0);
	auto time2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dtime = time2 - time;
	std::cerr << "RGB to Lab: " << dtime.count() << std::endl;
	time = std::chrono::high_resolution_clock::now();

	bilateral(im2, im2, sigma+0.5, sigma_col, sigma, W, H, mask);
	time2 = std::chrono::high_resolution_clock::now();
	dtime = time2 - time;
	std::cerr << "Bilateral filter: " << dtime.count() << std::endl;
	time = std::chrono::high_resolution_clock::now();

	#ifdef DEBUG_BIL
		for(int i = 0; i < W*H; i++) {
			if(mask == nullptr || mask[i]) {
				Color col = lab_to_rgb(im2[i]);
				for(int c = 0; c < 3; c++) im[3*i+c] = col.get(c);
			} else if(mask != nullptr) for(int c = 0; c < 3; c++) im[3*i+c] = 0;
		}
		stbi_write_png("bilateral.png", W, H, 3, im, 0);
		time2 = std::chrono::high_resolution_clock::now();
		dtime = time2 - time;
		std::cerr << "Lab to RGB and saving: " << dtime.count() << std::endl;
		time = std::chrono::high_resolution_clock::now();
	#endif

	return applySobelToBil(im2, W, H, mask, threshold, size_threshold);
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

void add_intersections(double r, double c, double s, int W, int H,
						const std::vector<std::pair<int, int>> &zone,
						std::vector<std::pair<int, int>> &ps) {
	int zs = zone.size();
	double hs = H*s, wc = W*c;
	ps.clear();
	if(zs < 4) {
		double hw = hs+wc;
		if(hs >= r) {
			ps.emplace_back(0, r/s);
			if(hw < r) ps.emplace_back((r - hs)/c, H);
		} else if(hw >= r) ps.emplace_back((r - hs)/c, H);
		if(wc >= r) {
			ps.emplace_back(r/c, 0);
			if(hw < r) ps.emplace_back(W, (r - wc)/s);
		} else if(hw >= r) ps.emplace_back(W, (r - wc)/s);
	} else {
		for(int i = 1; i < zs; i++) {
			double r0 = zone[i-1].first*c + zone[i-1].second*s;
			double r1 = zone[i].first*c + zone[i].second*s;
			double t = (r - r0) / (r1 - r0);
			if(0 <= t && t <= 1) {
				int x = zone[i].first*t + zone[i-1].first*(1-t);
				int y = zone[i].second*t + zone[i-1].second*(1-t);
				ps.emplace_back(x, y);
			}
		}
	}
	std::sort(ps.begin(), ps.end());
}

void saveHoughTransform(const char* filename, double *hough, int R, int T) {
	uchar *im = new uchar[R*T];
	double maxH = 0;
	for(int i = 0; i < R*T; i++)
		maxH = std::max(maxH, hough[i]);
	double f = 255.0 / maxH;
	for(int i = 0; i < R*T; i++)
		im[i] = hough[i] * f;
	stbi_write_png(filename, R, T, 1, im, 0);
}

PA::LinesData* PA::get_lines(PA::ProblemData* data, const std::vector<std::pair<int, int>> &zone) {
	// Time
	auto time = std::chrono::high_resolution_clock::now();

	// Important variables
	double W = data->W, H = data->H;
	double diag = std::sqrt(W*W + H*H);
	int R = 1.96 * std::pow(diag, 0.8);
	int T = 13.5 * std::pow(diag, 0.5);
	double r_step = diag / R;
	double t_step = 1.5 * M_PI / T;
	// Array of the Hough Transform
	double *res = new double[R*T];
	for(int i = 0; i < R*T; i++) res[i] = 0;

	// Compute the Hough Transform
	for(int x = 0; x < W; x++) {
		for(int y = 0; y < H; y++) {
			int pix = x + y*W;
			if(data->no[pix] < 0.001*data->m) continue;
			double ttm = std::atan2(y, x);
			double tt0 = (std::floor(ttm / t_step) + 0.5) * t_step - 0.5*M_PI;
			double tt1 = ttm + 0.5*M_PI;
			for(double tt = tt0; tt <= tt1; tt += t_step) {
				double diff_angle = std::abs(std::sin(tt - data->an[pix]));
				if(diff_angle > 0.8) continue;
				double c = std::cos(tt), s = std::sin(tt);
				double r = x*c + y*s;
				if(r < 0) continue;
				double val = (1.0 - std::pow(diff_angle, 0.66)) * (1.0 + 1.5 * data->no[pix]/data->m);
				int ri = int(r / r_step);
				double dr = (r - ri * r_step) / r_step;
				int ti = int((tt + M_PI/2.0) / t_step);
				if(dr > 0.5) {
					res[ri + ti * R] += (1.5 - dr) * val;
					if(ri+1 < R) res[ri+1 + ti * R] += (dr - 0.5) * val;
				} else {
					res[ri + ti * R] += (0.5 + dr) * val;
					if(ri-1 >= 0) res[ri-1 + ti * R] += (0.5 - dr) * val;
				}
			}
		}
	}
	#ifdef DEBUG_HOU
	saveHoughTransform("hough.png", res, R, T);
	#endif

	// Printing time
	auto time2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dtime = time2 - time;
	std::cerr << "Hough transform: " << dtime.count() << std::endl;
	time = std::chrono::high_resolution_clock::now();

	// Constante for the importance of the readjusment
	int x0, y0, x1, y1;
	if(zone.size() < 4) x0 = 0, y0 = 0, x1 = W, y1 = H;
	else {
		x0 = W, x1 = 0, y0 = H, y1 = 0;
		for(const auto p : zone)
			x0 = std::min(x0, p.first), x1 = std::max(x1, p.first), y0 = std::min(y0, p.second), y1 = std::max(y1, p.second);
	}
	double w2 = x1-x0, h2 = y1-y0;
	double norm = std::sqrt(w2*w2 + h2*h2);
	w2 /= norm, h2 /= norm;
	double cons = std::max(0.2, 1.25 - 1.4 * std::abs(w2 - h2) / double(w2 + h2));

	// Variables containing lines sorted by decreasing score
	std::vector<int> lines;
	for(int i = 0; i < R*T; i++) lines.push_back(i);
	std::sort(lines.begin(), lines.end(), [&res](int a, int b) { return res[a] > res[b]; });
	lines.resize(std::min(lines.size(), size_t(750 * std::pow(diag, 0.8))));
	int size_lines = lines.size();
	std::vector<PA::Line> ls;
	std::vector<int> indices;

	// Readjust score in function of the crossing area
	std::vector<double> cross_dis;
	std::vector<std::pair<int, int>> ps;
	for(int i : lines) {
		double r = (i % R + 0.5) * r_step, t = (int(i / R) + 0.5) * t_step - M_PI/2.0;
		double c = std::cos(t), s = std::sin(t);
		add_intersections(r, c, s, W, H, zone, ps);
		double dis = 0;
		int pss = ps.size();
		for(int i = 1; i < pss; i += 2) {
			double dx = ps[i].first - ps[i-1].first, dy = ps[i].second - ps[i-1].second;
			dis += std::sqrt(dx*dx + dy*dy);
		}
		cross_dis.push_back(dis);
	}
	double max_d = *std::max_element(cross_dis.begin(), cross_dis.end());
	for(int i = 0; i < size_lines; i++) res[lines[i]] *= 2.0 / (cons + cross_dis[i] / max_d);
	std::sort(lines.begin(), lines.end(), [&res](int a, int b) { return res[a] > res[b]; });

	// Removing and merging close lines
	int i = 0;
	double lim_dist = 0.0142 * diag;
	while(ls.size() < 200 && i < size_lines) {
		bool add = true;
		int ind = lines[i++];
		int ri = ind % R;
		int ti = ind / R;
		/*******************************/
		/* Smmothing with neighborhood */
		double min_neig = res[ind];
		double r_unw_sum = 0, t_unw_sum = 0;
		double r_w_sum = 0, t_w_sum = 0;
		double w = 0;
		double n = 0;
		for(int r2 = -1; r2 <= 1; r2 ++) {
			for(int t2 = -1; t2 <= 1; t2 ++) {
				int rb = ri+r2, tb = ti+t2;
				if(0 <= rb && rb < R && 0 <= tb && tb < T) {
					double val = res[rb + tb*R];
					min_neig = std::min(min_neig, val);
					r_unw_sum += r2, t_unw_sum += t2;
					r_w_sum += val*r2, t_w_sum += val*t2;
					w += val;
					n ++;
				}
			}
		}
		w -= n*min_neig;
		double r_add = (r_w_sum - r_unw_sum*min_neig) / w;
		double t_add = (t_w_sum - t_unw_sum*min_neig) / w;
		double r = (ri + r_add + 0.5) * r_step;
		double t = (ti + t_add + 0.5) * t_step - M_PI/2;
		/*******************************/
		double c = std::cos(t), s = std::sin(t);
		add_intersections(r, c, s, W, H, zone, ps);
		double x0 = ps[0].first, y0 = ps[0].second;
		double x1 = ps.back().first, y1 = ps.back().second;
		double lim_dist2 = (0.35 + std::abs(c*w2) + std::abs(s*h2)) / 1.35 * lim_dist;
		std::vector<int>::iterator ind_it = indices.begin();
		for(PA::Line &l : ls) {
			int dth = std::abs(t - l.theta);
			if((dth > 0.42) && !(r < lim_dist && l.rho < lim_dist && std::abs(dth - M_PI) > 0.42)) continue;
			double co = std::cos(l.theta), si = std::sin(l.theta);
			double d0 = co * x0 + si * y0 - l.rho;
			double d1 = co * x1 + si * y1 - l.rho;
			bool cross = (d0 < 0) != (d1 < 0);
			double dm = cross ? (d0*d0 + d1*d1) / 2.0 / std::abs(d0 - d1) : std::abs(d0 + d1) / 2.0;
			if(dm < lim_dist2) {
				add = false;
				if(dm > 0.25*lim_dist2) break;
				double alpha = 10e4*cons/std::pow(diag, 1.8) * (1.0 - dm / (0.25*lim_dist2)) * res[ind] / (res[ind] + res[*ind_it]);
				res[*ind_it] += res[ind];
				l.set((1 - alpha) * l.rho + alpha * r, (1 - alpha) * l.theta + alpha * t, true);
				break;
			}
			ind_it ++;
		}
		if(!add) continue;
		ls.emplace_back(r, t, true);
		indices.push_back(ind);

		// double max_d = 0.008 * diag;
		// double co = std::cos(t), si = std::sin(t);
		// std::vector<int> ps;
		// Line l(r, t, true);
		// for(int x = 0; x < data->W; x++) {
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
		// 		if(data->no[pix] > 0 && std::abs(std::sin(t - data->an[pix])) < 0.25)
		// 			ps.push_back(pix);
		// 		if(y >= y0) y ++;
		// 		else y --;
		// 	}
		// }
		// Line l2 = pca(ps, data->no, data->an, data->W);
		// ls.push_back(l2);
	}
	time2 = std::chrono::high_resolution_clock::now();
	dtime = time2 - time;
	std::cerr << "Finding lines: " << dtime.count() << std::endl;

	return new PA::LinesData(R, T, res, ls);
}
