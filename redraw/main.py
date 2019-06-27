import argparse
import sys
import pylab as pl
import numpy as np

class Line:
	def __init__(self, x1, y1, x2, y2, col=(0, 0, 0)):
		self._x = [x1, x2]
		self._y = [y1, y2]
		self._col = col
		self._a = (y2 - y1) / (x2 - x1) if x2 != x1 else float('inf')
		self._b = y1 - self._a * x1
		self._theta = - np.math.atan(1.0 / self._a)
		self._c = np.cos(self._theta)
		self._s = np.sin(self._theta)
		self._r = self._c * x1 + self._s * y1
	
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
		
	def plot(self):
		pl.plot(self._x, self._y, color=self._col)

# Will be probably removed
def to_vanish(li, c):
	d1 = abs(c[0] - li._x[0]) + abs(c[1] - li._y[0])
	d2 = abs(c[0] - li._x[1]) + abs(c[1] - li._y[1])
	i = 0 if d1 < d2 else 1
	a = (c[1] - li._y[i]) / (c[0] - li._x[i])
	if abs(li._x[1-i]) > abs(li._y[1-i]):
		x, y = li._x[1-i], c[1] + a * (li._x[1-i] - c[0])
	else:
		x, y = c[0] + (li._y[1-i] - c[1]) / a, li._y[1-i]
	return Line(x, y, li._x[i], li._y[i], li._col)

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

	# Creation of Figure 2
	pl.figure(2)
	pl.xlim(0, W)
	pl.ylim(H, 0)

	for g in groups:
		if g not in dots:
			ma = sum(l._a for l in groups[g]) / len(groups[g])
			mb = sum(l._b for l in groups[g]) / len(groups[g])
			print(ma, mb)
			pl.plot([0, W], [mb, mb+W*ma])
			vp = horizon.inter(Line(0, mb, W, mb+W*ma))
			print(vp)
		else:
			c = dots[g]
			x, y = horizon.proj(*c)
			xls = [l.x(im_y+im_h) for l in groups[g]]
			x0, x1 = min(xls), max(xls)
			nls = len(xls)
			col = groups[g][0]._col
			for i in range(nls):
				l = Line(x, y, x0 + (x1 - x0) * i / (nls-1), im_y+im_h, col)
				l.plot()
	pl.show()

if __name__ == "__main__":
	main()