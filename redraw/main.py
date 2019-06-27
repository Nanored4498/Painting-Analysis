import argparse
import sys
import pylab as pl
import numpy as np

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

	# Recompute vanishing points
	gs = sorted(dots.keys(), key=lambda g: dots[g][0])
	dx = 0.5 * (dots[gs[2]][0] - dots[gs[0]][0])
	a = (dots[gs[0]][0] + dots[gs[1]][0] + dots[gs[2]][0]) / 3.0 - dx
	dy = 0.5 * (dots[gs[2]][1] - dots[gs[0]][1])
	b = (dots[gs[0]][1] + dots[gs[1]][1] + dots[gs[2]][1]) / 3.0 - dy
	dots2 = {gs[i] : (a+i*dx, b+i*dy) for i in range(3)}
	horizon2 = Line(a, b, a+dx, b+dy)
	bottom = horizon2.translate(0, H - 0.5 * (horizon2.y(im_x) + horizon2.y(im_x+im_w)))

	# Recompute parallel lines
	datas = {}
	dx, dy = 0, 0
	for g in dots2:
		x, y = dots2[g]
		im_bot = im_y+im_h
		xls = [l.x(im_bot) for l in groups[g]]
		x0, x1 = min(xls), max(xls)
		l0, l1 = Line(x, y, x0, im_bot), Line(x, y, x1, im_bot)
		p0, p1 = bottom.inter(l0), bottom.inter(l1)
		nls = len(xls)
		col = groups[g][0]._col
		a, b = p0
		dx += (p1[0] - p0[0]) / (nls-1) / 3.0
		dy += (p1[1] - p0[1]) / (nls-1) / 3.0
		datas[g] = x, y, a, b, col, nls
	d0c, d0s = 0, 0
	for g in datas:
		a = datas[g][2]
		da = (a - int(a / dx) * dx) / dx * 2 * np.pi
		d0c += np.math.cos(da)
		d0s += np.math.sin(da)
	t0 = np.math.atan2(d0s, d0c)
	if t0 < 0 : t0 += 2 * np.pi
	d0 = t0 * dx / (2 * np.pi)
	for g in datas:
		x, y, a, b, col, nls = datas[g]
		da = a - int(a / dx) * dx
		if 2*(da - d0) > dx: da -= dx
		elif 2*(d0 - da) > dx: da += dx
		a += d0 - da
		b += (d0 - da) * dy / dx
		for i in range(nls):
			l = Line(x, y, a + i * dx, b + i * dy, col)
			l.plot()
	pl.show()

if __name__ == "__main__":
	main()