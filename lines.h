#ifndef LINES_H
#define LINES_H

#include <cmath>
#include <vector>

#define MODE_PI
#define NUM_PASS 6

#define uchar unsigned char

namespace PA {

struct Line {
	double a, b, rho, theta;

	Line(double x, double y, bool hough=false) {
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
	double m;

	ProblemData(int W, int H, double *no, double *an, double m):
		W(W), H(H), no(no), an(an), m(m) {}

	~ProblemData() {
		delete [] no;
		delete [] an;
	}
};

ProblemData* applySobel(uchar* im, int W, int H);
void save_sobel(const char *filename, PA::ProblemData *data);
std::vector<Line> get_lines(ProblemData *data);

}

#undef uchar

#endif // LINES_H
