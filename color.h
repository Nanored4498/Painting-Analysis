#include <cmath>
#include <iostream>

#define Xn 95.0
#define Yn 100.0
#define Zn 108.9

#define Delta (6.0 / 29.0)
const static double delta2 = std::pow(Delta, 2.0);
const static double delta3 = std::pow(Delta, 3.0);

struct Color {
	double cs[3];

	Color(double x=0, double y=0, double z=0) {
		cs[0] = x;
		cs[1] = y;
		cs[2] = z;
	}

	inline double x() const { return cs[0]; }
	inline double y() const { return cs[1]; }
	inline double z() const { return cs[2]; }
	inline double get(int i) { return cs[i]; }

	Color operator+(const Color &other) const { return Color(x()+other.x(), y()+other.y(), z()+other.z()); }
	Color operator-(const Color &other) const { return Color(x()-other.x(), y()-other.y(), z()-other.z()); }
	Color& operator+=(const Color &other) { return (*this = *this + other); }
	Color operator*(double d) const { return Color(x()*d, y()*d, z()*d); }
	Color operator/(double d) const { return Color(x()/d, y()/d, z()/d); }
};

inline Color operator*(double d, const Color &c) { return c * d; }

inline std::ostream& operator<<(std::ostream &stream, const Color &c) {
	return stream << "(" << c.x() << ", " << c.y() << ", " << c.z() << ")";
}

inline double f_xyz_lab(double t) {
	if(t > delta3) return std::pow(t, 1.0/3.0);
	return t / (3 * delta2) + 4.0 / 29.0;
}

inline Color xyz_to_lab(const Color &c) {
	double L = 116 * f_xyz_lab(c.y() / Yn) - 16;
	double a = 500 * (f_xyz_lab(c.x() / Xn) - f_xyz_lab(c.y() / Yn));
	double b = 200 * (f_xyz_lab(c.y() / Yn) - f_xyz_lab(c.z() / Zn));
	return Color(L, a, b);
}

inline Color rgb2_to_xyz(const Color &c) {
	double X = 0.4124564 * c.x() + 0.3575761 * c.y() + 0.1804375 * c.z();
	double Y = 0.2126729 * c.x() + 0.7151522 * c.y() + 0.0721750 * c.z();
	double Z = 0.0193339 * c.x() + 0.1191920 * c.y() + 0.9503041 * c.z();
	return Color(X, Y, Z) * 100.0;
}

inline Color rgb_to_rgb2(const Color &c) {
	double r = c.x() / 255.0;
	double g = c.y() / 255.0;
	double b = c.z() / 255.0;
	double R = r <= 0.04045 ? r / 12.92 : std::pow((r + 0.055) / 1.055, 2.4);
	double G = g <= 0.04045 ? g / 12.92 : std::pow((g + 0.055) / 1.055, 2.4);
	double B = b <= 0.04045 ? b / 12.92 : std::pow((b + 0.055) / 1.055, 2.4);
	return Color(R, G, B);
}

inline Color rgb_to_xyz(const Color &c) {
	return rgb2_to_xyz(rgb_to_rgb2(c));
}

inline Color rgb_to_lab(const Color &c) {
	return xyz_to_lab(rgb_to_xyz(c));
}

inline Color rgb2_to_rgb(const Color &c) {
	double R = c.x() <= 0.0031308 ? c.x() * 12.92 : std::pow(c.x(), 1.0 / 2.4) * 1.055 - 0.055;
	double G = c.y() <= 0.0031308 ? c.y() * 12.92 : std::pow(c.y(), 1.0 / 2.4) * 1.055 - 0.055;
	double B = c.z() <= 0.0031308 ? c.z() * 12.92 : std::pow(c.z(), 1.0 / 2.4) * 1.055 - 0.055;
	return Color(R, G, B) * 255.0;
}

inline Color xyz_to_rgb2(const Color &c) {
	double x = c.x() / 100.0;
	double y = c.y() / 100.0;
	double z = c.z() / 100.0;
	double R =   3.2404542 * x - 1.5371385 * y - 0.4985314 * z;
	double G = - 0.9692660 * x + 1.8760108 * y + 0.0415560 * z;
	double B =   0.0556434 * x - 0.2040259 * y + 1.0572252 * z;
	return Color(R, G, B);
}

inline double f_lab_xyz(double t) {
	if(t > Delta) return std::pow(t, 3.0);
	return (t - 4.0 / 29.0) * 3 * delta2;
}

inline Color lab_to_xyz(const Color &c) {
	double l = (c.x() + 16) / 116.0;
	double Y = Yn * f_lab_xyz(l);
	double X = Xn * f_lab_xyz(c.y() / 500.0 + l);
	double Z = Zn * f_lab_xyz(-c.z() / 200.0 + l);
	return Color(X, Y, Z);
}

inline Color lab_to_rgb2(const Color &c) {
	return xyz_to_rgb2(lab_to_xyz(c));
}

inline Color lab_to_rgb(const Color &c) {
	return rgb2_to_rgb(lab_to_rgb2(c));
}