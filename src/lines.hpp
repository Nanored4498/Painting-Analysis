#ifndef LINES_H
#define LINES_H

#include <cmath>
#include <vector>
#include "color.hpp"

#define MIN_SOBEL_INTENSITY 75

#define MODE_PI

#define uchar unsigned char

namespace PA {

struct Line {
	double a, b, rho, theta;

	Line(double x, double y, bool hough=false) {
		set(x, y, hough);
	}

	void set(double x, double y, bool hough=false) {
		if(hough) {
			rho = x;
			theta = y;
			a = - 1.0 / std::tan(theta);
			b = rho / std::sin(theta);
		} else {
			a = x;
			b = y;
			double t = - std::atan(1.0 / a);
			theta = b > 0 ? (t >= 0 ? t : t + M_PI) : (t <= 0 ? t : t - M_PI);
			rho = b * std::sin(theta);
		}
	}
};

struct ProblemData {
	int W, H;
	double *no, *an;
	Color *bil;
	double m;

	ProblemData(int W, int H, double *no, double *an, Color *bil, double m):
		W(W), H(H), no(no), an(an), bil(bil), m(m) {}

	~ProblemData() {
		delete[] no;
		delete[] an;
		delete[] bil;
	}
};

struct LinesData {
	int R, T;
	double *score;
	std::vector<Line> lines;

	LinesData(int R, int T, double *score, std::vector<Line> lines):
		R(R), T(T), score(score), lines(lines) {}

	~LinesData() {
		delete[] score;
	}
};

ProblemData* applySobel(uchar* im, int W, int H,
						bool *mask=nullptr,
						double threshold=6.7,
						double size_threshold=7.5);
ProblemData* applySobelToBil(Color* im, int W, int H,
						bool *mask=nullptr,
						double threshold=6.7,
						double size_threshold=7.5);

uint pixelColor(int x, int y, PA::ProblemData *data);
void save_sobel(const char *filename, PA::ProblemData *data);

LinesData* get_lines(ProblemData *data, const std::vector<std::pair<int, int>> &zone);

}

#undef uchar

#endif // LINES_H
