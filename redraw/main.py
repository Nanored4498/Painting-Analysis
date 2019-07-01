import argparse
import sys
import pylab as pl
import numpy as np
import sympy

class Line:
	def __init__(self, x1, y1, x2, y2, col=(0, 0, 0)):
		self._x = np.array([x1, x2])
		self._y = np.array([y1, y2])
		self._col = col
		self._a = (y2 - y1) / (x2 - x1) if x2 != x1 else float('inf')
		self._b = y1 - self._a * x1
		self._compute_hough()
	
	def _compute_hough(self):
		self._theta = - np.math.atan(1.0 / self._a)
		self._c = np.cos(self._theta)
		self._s = np.sin(self._theta)
		self._r = self._c * self._x[0] + self._s * self._y[0]

	def dist(self, x, y):
		return abs(self._c * x + self._s * y - self._r)
	
	def proj(self, x, y):
		d = self._c * x + self._s * y - self._r
		return x - d * self._c, y - d * self._s
	
	def y(self, x): return self._b + self._a * x
	def x(self, y): return (y - self._b) / self._a
	
	def inter(self, other):
		x = (self._b - other._b) / (other._a - self._a)
		return (x, self.y(x))
	
	def translate(self, dx, dy):
		return Line(self._x[0]+dx, self._y[0]+dy, self._x[1]+dx, self._y[1]+dy)
		
	def plot(self):
		pl.plot(self._x, self._y, color=self._col)

def read_attributes(l, A):
	"""
	Read the all the attributes of `A` in the line `l` and return a dictionary
	containing the value of each attribute
	"""
	res = {}
	for a in A:
		b = a + '="'
		start = l.find(b) + len(b)
		end = l.find('"', start)
		res[a] = float(l[start:end])
	return res

def main():
	parser = argparse.ArgumentParser(description='Redraw the perspective lines of a SVG correctly')
	parser.add_argument('filename', type=str,
						help='The path to the SVG file to redraw')
	args = parser.parse_args()

	# Opening the file
	try:
		f = open(args.filename, "r")
	except FileNotFoundError:
		print("File not found:", args.filename, file=sys.stderr)
		exit(1)

	# Reading the size
	l = f.readline().replace(' ', '')
	size = read_attributes(l, ['width', 'height'])
	W, H = size['width'], size['height']
	
	# Reading painting position and size
	l = f.readline().replace(' ', '')
	a = read_attributes(l, ['x', 'y', 'width', 'height'])
	im_x, im_y = a['x'], a['y']
	im_w, im_h = a['width'], a['height']
	im_bot = im_y+im_h
	im_right = im_x+im_w

	# Creation of Figure 1
	pl.figure(1)
	pl.xlim(0, W)
	pl.ylim(H, 0)
	
	# Reading groups of parallel lines plus the horizon
	groups = {}
	while True:
		l = f.readline().replace(' ', '').replace('\t', '')
		if not l.startswith("<line"): break
		cs = read_attributes(l, ['x1', 'x2', 'y1', 'y2'])
		start = l.find('rgb(') + 4
		end = l.find(',', start)
		r = int(l[start:end]) / 255.0
		start = end+1
		end = l.find(',', start)
		g = int(l[start:end]) / 255.0
		start = end+1
		end = l.find(')', start)
		b = int(l[start:end]) / 255.0
		li = Line(cs['x1'], cs['y1'], cs['x2'], cs['y2'], (r, g, b))
		g = int(255 * (r + 256*(g + 256*b)))
		if g in groups: groups[g].append(li)
		else: groups[g] = [li]
		li.plot()
	horizon = groups[g][0]
	del groups[g]
	
	# Reading dots and making the correspondance whith groups of lines
	dots = {}
	while True:
		if not l.startswith("<circle"): break
		cs = read_attributes(l, ['cx', 'cy'])
		best = float('inf')
		bg = 0
		x, y = cs['cx'], cs['cy']
		for g in groups:
			li = groups[g][0]
			dis = li.dist(x, y)
			if dis < best:
				best = dis
				bg = g
		dots[bg] = (x, y)
		pl.scatter(x, y)
		l = f.readline().replace(' ', '').replace('\t', '')

	# Need 3 vanishing points
	assert(len(dots) == 3)

	# Creation of Figure 2
	pl.figure(2)
	pl.xlim(0, W)
	pl.ylim(H, 0)

	################################################
	# Importance of the vanish points of depth lines
	center_coeff = 5
	################################################
	
	# Computing the mean of horizontal lines coefficient
	for g in groups:
		if g not in dots:
			gh = g
			ma = sum(l._a for l in groups[g]) / len(groups[g])

	# Compute a new horizon and a bottom line
	gs = sorted(dots.keys(), key=lambda g: dots[g][0])
	dx = 0.5 * (dots[gs[2]][0] - dots[gs[0]][0])
	a = (dots[gs[0]][0] + center_coeff * dots[gs[1]][0] + dots[gs[2]][0]) / (2 + center_coeff) - dx
	dy = 0.5 * (dots[gs[2]][1] - dots[gs[0]][1])
	b = (dots[gs[0]][1] + center_coeff * dots[gs[1]][1] + dots[gs[2]][1]) / (2 + center_coeff) - dy
	horizon2 = Line(a-3*dx, b-3*dy, a+3*dx, b+3*dy, horizon._col)
	y_translation = im_bot - 0.5 * (horizon2.y(im_x) + horizon2.y(im_right))
	d_translation = dx * y_translation / (dx**2 + dy**2) ** 0.5
	bottom = horizon2.translate(0, y_translation)

	# Initial values to optimize
	xs = [dots[gs[i]][0] for i in range(3)]
	ys = [dots[gs[i]][1] for i in range(3)]
	ds = []
	# Coefficients in the objective functions
	cxs = [(xs[1]-xs[0])**(-2), center_coeff * 4*(xs[2]-xs[0])**(-2), (xs[2]-xs[1])**(-2)]
	cys = [(ys[1]-ys[0])**(-2), center_coeff * 4*(ys[2]-ys[0])**(-2), (ys[2]-ys[1])**(-2)]
	p1s, p2s = [], []
	nls = []

	# Computing some usefull values
	for g in gs:
		x, y = dots[g]
		ols = sorted(groups[g], key=lambda l: l.x(im_bot))
		l0, l1 = Line(x, y, ols[0].x(im_bot), im_bot), Line(x, y, ols[-1].x(im_bot), im_bot)
		p0, p1 = bottom.inter(l0), bottom.inter(l1)
		nl = len(ols)
		col = groups[g][0]._col
		dx = (p1[0] - p0[0]) / (nl-1)
		dy = (p1[1] - p0[1]) / (nl-1)
		nls.append(nl)
		ds.append((dx**2 + dy**2) ** 0.5)
		for l, ps in [(ols[0], p1s), (ols[-1], p2s)]:
			x = l.x(im_bot)
			if x < im_x: x, y = im_x, l.y(im_x)
			elif x > im_right: x, y = im_right, l.y(im_right)
			else: y = l.y(x)
			ps.append((x, y))
	cds = (min(ds) / 3) ** (-2)

	# V = sympy.symbols('x0 x2 y0 y2 d0 d2')
	# X0, X2, Y0, Y2, D0, D2 = V
	# X1 = (X0*D2 + X2*D0) / (D0 + D2)
	# Y1 = (Y0*D2 + Y2*D0) / (D0 + D2)
	# D1 = 2 * D0 * D2 / (D0 + D2)
	# NU = sympy.sqrt((X2 - X1)**2 + (Y2 - Y1)**2)
	# U = (X2 - X1) / NU, (Y2 - Y1) / NU
	# N = (Y1 - Y2) / NU, (X2 - X1) / NU
	# X, Y, D = [X0, X1, X2], [Y0, Y1, Y2], [D0, D1, D2]
	# BS, CS = [], []
	# F = 0
	# for i in range(3):
	# 	APX, APY = p1s[i][0] - X[i], p1s[i][1] - Y[i]
	# 	N_AP = sympy.simplify(N[0] * APX + N[1] * APY)
	# 	BX = sympy.simplify(X[i] + APX * d_translation / N_AP)
	# 	BY = sympy.simplify(Y[i] + APY * d_translation / N_AP)
	# 	CX = BX + (nls[i] - 1) * D[i] * U[0]
	# 	CY = BY + (nls[i] - 1) * D[i] * U[1]
	# 	BS.append((BX, BY))
	# 	CS.append((CX, CY))
	# 	ACX, ACY = CX - X[i], CY - Y[i]
	# 	CP2X, CP2Y = p2s[i][0] - CX, p2s[i][1] - CY
	# 	DP2 = (CP2X * ACY - CP2Y * ACX) ** 2 / (ACX**2 + ACY**2)
	# 	F += cxs[i]*(X[i]-xs[i])**2 + cys[i]*(Y[i]-ys[i])**2 + cds * DP2
	# G = [F.diff(v) for v in V]
	# H = [[g.diff(v) for v in V] for g in G]

	# def newton(v2):
	# 	s = [(V[i], v2[i]) for i in range(6)]
	# 	fv = float(F.subs(s))
	# 	gv = np.array([float(g.subs(s)) for g in G])
	# 	print(fv, max(abs(gv)))
	# 	hv = np.array([[float(hi.subs(s)) for hi in h] for h in H])
	# 	hvi = np.linalg.inv(hv)
	# 	return v2 - hvi.dot(gv)
	
	# v2 = np.array([xs[0], xs[2], ys[0], ys[2], ds[0], ds[2]])
	# s = [(V[i], v2[i]) for i in range(6)]
	# xs, ys, ds = [float(x.subs(s)) for x in X], [float(y.subs(s)) for y in Y], [float(d.subs(s)) for d in D]
	# print(xs)
	# print(ys)
	# print(ds)
	# for i in range(3):
	# 	v2 = newton(v2)
	# 	s = [(V[i], v2[i]) for i in range(6)]
	# 	xs, ys, ds = [float(x.subs(s)) for x in X], [float(y.subs(s)) for y in Y], [float(d.subs(s)) for d in D]
	# 	print(xs)
	# 	print(ys)
	# 	print(ds)
	
	xs = [136.64530420611752, 1150.345165488534, 2522.897642552893]
	ys = [819.5731053628471, 764.5395970563529, 690.0240712753191]
	ds = [99.95955861425124, 114.99181170194613, 135.34552486638526]

	# Recompute the new horizon and the bottom line
	dx, dy = xs[2] - xs[0], ys[2] - ys[0]
	horizon2 = Line(xs[0]-dx, ys[0]-dx, xs[2]+dx, ys[2]+dx, horizon._col)
	nd = (dx**2 + dy**2) ** 0.5
	dx /= nd
	dy /= nd
	print(dx, dy)
	y_translation = d_translation / dx
	bottom = horizon2.translate(0, y_translation)

	for i in range(3):
		g = gs[i]
		col = groups[g][0]._col
		x, y = xs[i], ys[i]
		px, py = p1s[i]
		apx, apy = px - x, py - y
		t = d_translation / (dx*apy - dy*apx)
		px, py = x + apx*t, y + apy*t
		d = ds[i]
		for _ in range(nls[i]):
			l = Line(x, y, px, py, col)
			l.plot()
			px += dx * d
			py += dy * d

	# for l in groups2[gs[2]]:
	# 	x, y = l.inter(groups2[gs[0]][0])
	# 	l2 = Line(x-20*dx, y-20*dy, x+20*dx, y+20*dy, groups[gh][0]._col)
	# 	l2.plot()

	# Show figures
	pl.show()

if __name__ == "__main__":
	main()